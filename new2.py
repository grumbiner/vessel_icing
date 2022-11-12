import sys
import math
import numpy as np
import numpy.ma as ma

import pygrib

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg') #for batch mode
#matplotlib.use('Qt5Agg') #for interactive mode
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator

import cartopy.crs as ccrs
import cartopy.feature as cfeature

#pythonpath = mmablib/py
from struct import *
from const import *
from ijpt  import *
from latpt import *
from grid  import *

#---------------------------------------------------------------
#Vessel Icing constants, cm/hr:
icing_a = 2.73E-02
icing_b = 2.91E-04
icing_c = 1.84E-06
#freezing point of sea water in algorithm
tf = -1.7

#working variables:
mapper  = global_nmc_nthdeg(4) 
nx = int(mapper.nx)
ny = int(mapper.ny)

#debug print("nx etc. ",nx, ny)
#exit(0)

#read static sst, land, ice
#RG: args:
res="0p25"
tag="20221111"

#RG: Could also manage grib gets from urllib3 

hr="000"

fname="sst."+res+"."+tag+".f"+hr+".grib2"
grbsst = pygrib.open(fname)
#debug print("grbsst = ",grbsst, flush=True)

for x in grbsst:
  #debug print(x.shortName, x.name, x.level, x.typeOfLevel, flush=True)
  sst = x.values
  print(x.shortName, sst.max(), sst.min(), flush=True)

#Get lats, lons from the grib message itself:
lats,lons = x.latlons()

grbsst.close()

#convert to C if needed
if (sst.max() > 150) :
   sst -= 273.15

#masks based on sst:
om3 = ma.masked_array(sst > 12.0)
om4 = ma.masked_array(sst < tf - 1.)
omask = ma.mask_or(om3, om4)

#exit(0)

#------- ice and land:
fname = "landice."+res+"."+tag+".f"+hr+".grib2"
grblandice = pygrib.open(fname)
#debug print("landice = ",grblandice, flush=True)

for x in grblandice:
  #debug print(x.shortName, x.name, x.level, x.typeOfLevel, flush=True)
  y = x.values
  if (x.shortName == "lsm"):
    land = y
  elif (x.shortName == "ci"):
    ice  = y
  print(x.shortName, y.max(), y.min(), flush=True)

grblandice.close()

#exit(0)

#----------------------------------------------------
#  Manage masks
#OMask = true if it not possible to have icing:
#omask = ma.masked_array(land > 0.5 or ice > 0.5 or sst < tf or sst > 12.0)
om1 = ma.masked_array(land > 0.5)
om2 = ma.masked_array(ice > 0.5)
omask = ma.mask_or(omask, om1)
omask = ma.mask_or(omask, om2)

#mask = true if it is possible to have icing -- negate omask
mask = ~omask
i_index = mask.nonzero()[1]
j_index = mask.nonzero()[0]
npts = len(i_index)
#debug print("npts = ",npts, flush=True)

#----------------------------------------------------
# prepare variables
ll = latpt
icing_rate  = np.zeros((ny, nx))

histogram = np.zeros(8000)
rate_axis = np.zeros(len(histogram))

sum = 0.0
sumsq = 0.0
sumarea = 0.0
icing_plus = 0.0

dtime = int(3)
#-------------------------------------------------------------------
for hour in range(0, 80*dtime + 1, dtime):
#for hour in range(0, 8*dtime + 1, dtime):
  hr = "{:03d}".format(hour)

  fname="running."+res+"."+tag+".f"+hr+".00.grib2"
  grbs = pygrib.open(fname) 
  #debug print(hour, fname, grbs, flush=True)

  for x in grbs:
    #debug print(x.shortName, x.name, x.level, x.typeOfLevel, flush=True)
    y = x.values
    if (x.shortName == "2t"):
      t2m = y
    elif (x.shortName == "10u"):
      u10 = y
    elif (x.shortName == "10v"):
      v10 = y
    else:
      print("unknown/unused variable ",x.shortName, flush=True)
    #debug: print(x.shortName, y.max(), y.min(), flush=True)

# transform units, derive quantities:
  if (t2m.max() > 150):
    t2m -= 273.15
  v10 *= v10
  u10 *= u10
  speed = (u10+v10)
#  speed = np.sqrt(speed)

# main loop:
  for k in range (0,npts):
      i = i_index[k]
      j = j_index[k]

      if (t2m[j,i] > 0 or t2m[j,i] < -40.0 ) :
        icing_rate[j,i] = 0.0
        icing_plus = 0.0
      else :
        speed[j,i] = sqrt(speed[j,i])
        PR = speed[j,i]*(tf - t2m[j,i])/(1.+.4*(sst[j,i] - tf))
        icing_rate[j,i] = PR*(icing_a + icing_b*PR + icing_c*PR*PR)
        PRplus = (2.+speed[j,i])*(tf - (t2m[j,i]-2.) )/(1.+.4*( max(tf, (sst[j,i]-0.5)) - tf))
        icing_plus = PRplus*(icing_a + icing_b*PRplus + icing_c*PRplus*PRplus)
        #dI_dP = icing_a + 2.*icing_b*PR + 3.*icing_c*PR*PR
        #dPdw = (1./speed[j,i])
        #dPdTa = -1./(tf-t2m[j,i])
        #dPdTo =  -0.4/(1+0.4*(sst[j,i] - tf)) 
# cm/hr
      if (icing_plus > 0 ) :
        #print ("grid ", "{:3d}".format(hour), "{:7.3f}".format(ll.lat), 
        #          "{:7.3f}".format(ll.lon), "{:6.2f}".format(icing_rate[j,i]), 
        #          "{:6.2f}".format(icing_plus) )

        if (icing_rate[j,i] >= 0):
          irate = round(icing_rate[j,i]*10.0)
          da    = mapper.cellarea(j)
          sum     += icing_rate[j,i] * da
          sumsq   += icing_rate[j,i] * icing_rate[j,i] * da
          sumarea += da
          histogram[int(irate)] += da

#From pcolormesh sample ------------------------------------------
  levels = (0,0.01,0.7,2.0, 4.0, 7.5, 13.5, len(histogram)/10.)
#  cmap = plt.get_cmap('PiYG')
# hand craft a color map
  none = np.array([1., 1., 1., 1.])
  trace = np.array([0.2, 0.2, 0.2, 1.])
  light = np.array([0.5, 0.5, 1., 1])
  medium  = np.array([.5, 0., 0., 1])
  heavy   = np.array([1., 0., 0., 1])
  extreme = np.array([0., 0., 1., 1])
  vheavy  = np.array([1., 1., 0., 1])
  vvheavy = np.array([0., 1., 1., 1])
  cmap = ListedColormap([none,trace, light, medium, heavy, extreme, vheavy, vvheavy])

  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

  #proj = ccrs.LambertConformal(central_longitude = 150., central_latitude = 60., cutoff=25.)
  proj = ccrs.PlateCarree()
  ax  = plt.axes(projection = proj)
  fig = plt.figure(figsize=(12,9))
  ax  = fig.add_subplot(1,1,1, projection = proj)

#cartopy is not good about crossing 180 or 0

  #Barents: ax.set_extent((   0,  90, 50, 90), crs = ccrs.PlateCarree() )
  #Okhotsk: ax.set_extent(( 115, 165, 35, 65), crs = ccrs.PlateCarree() )
  #Bering:  ax.set_extent((-180,-130, 45, 90), crs = ccrs.PlateCarree() )
  #N. Atl:  ax.set_extent(( -95, -45, 40, 80), crs = ccrs.PlateCarree() )
  #Nordic: 
  ax.set_extent((-45,-1, 40, 90), crs = ccrs.PlateCarree() )

  #ax.gridlines(crs=ccrs.PlateCarree(),
  #             #xlocs=[180, 190, 200, 220, 240],
  #             xlocs=[150, 160, 170],
  #             ylocs=[40,50,60,66.6,70,80] )
  ax.gridlines(crs=ccrs.PlateCarree() )

  #ax.add_feature(cfeature.GSHHSFeature(levels=[1,2,3,4], scale="l") )
  ax.add_feature(cfeature.GSHHSFeature(levels=[1,2], scale="l") )
  #ax.coastlines()

  #cs = ax.pcolormesh(lons, lats, icing_rate, cmap=cmap, norm=norm, alpha = .990)
  cs = ax.pcolormesh(lons, lats, icing_rate, cmap=cmap, norm=norm)
  cb = plt.colorbar(cs)
  cbarlabel = '%s' % ("icing rates (cm/hr)")
  cb.set_label(cbarlabel, fontsize = 12)
  
  ax.set(xlabel = "Longitude")
  ax.set(ylabel = "Latitude")
  ax.set_title("Icing Rate (cm/hr) -- " + "{:03d}".format(hour) + " hour lead forecast from "+tag)

  #for interactive: plt.show()
  plt.savefig("b"+"{:03d}".format(hour) +".png")
  plt.close()

#----------------------------------------------------------------------

print ("average, rms of nonzero = ",sum / sumarea, sqrt(sumsq/sumarea))
histogram /= sumarea
light    = 0.0 # to 0.7 cm/hr
moderate = 0.0 # to 2.0 cm/hr
heavy    = 0.0 # to 4.0 cm/hr
extreme  = 0.0 # over 4.0 cm/hr (std.)
vheavy   = 0.0
vvheavy  = 0.0

for i in range (0,len(histogram) ):
  rate = float(i) / 10.0
  rate_axis[i] = rate
  if (rate < 0.7):
    light += histogram[i]
  elif (rate < 2.0):
    moderate += histogram[i]
  elif (rate < 4.0):
    heavy += histogram[i]
  elif (rate < 7.5):
    vheavy += histogram[i]
  elif (rate < 13.5):
    vvheavy += histogram[i]
  else:
    extreme += histogram[i]

  if (histogram[i] > 0):
    print ("hist", i/10.0, histogram[i])

print ("low-heavy ", light, moderate, heavy, light+moderate+heavy)
print ("vh-extreme ", vheavy, vvheavy, extreme)
print ("all ", light + moderate + heavy + vheavy + vvheavy + extreme)

#----------------------------------------------------------------------
fig, ax = plt.subplots()
ax.set(xlabel = "Icing rate (cm/hr)")
ax.set(ylabel = "Fraction of area with any icing")
ax.plot(rate_axis[0:135], histogram[0:135])
ax.grid()
plt.savefig("hist.png")
plt.close()
