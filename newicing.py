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
from matplotlib.ticker import MaxNLocator

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
#mapper  = global_nmc_nthdeg(1) 
mapper  = global_nmc_nthdeg(4) 
nx = int(mapper.nx)
ny = int(mapper.ny)
dx = mapper.dlon
dy = mapper.dlat
lons = np.mgrid[slice(dx/2, 360., dx)]
lats = np.mgrid[slice(90., -90.1, dy)]
print("nx etc. ",nx, ny, dx, dy)

ll = latpt
icing_rate  = np.zeros((ny, nx))

#read static sst, land, ice
res="0p25"
hr="000"
tag="20221111"

fname="sst."+res+"."+tag+".f"+hr+".grib2"
grbsst = pygrib.open(fname)
print("grbsst = ",grbsst)

for x in grbsst:
  #debug print(x.shortName, x.name, x.level, x.typeOfLevel, flush=True)
  sst = x.values
  print(sst.max(), sst.min())

#pygrib.close(fname)
grbsst.close()


#convert to C if needed
if (sst.max() > 150) :
   sst -= 273.15

#masks based on sst:
om3 = ma.masked_array(sst > 12.0)
om4 = ma.masked_array(sst < tf)
omask = ma.mask_or(om3, om4)

#exit(0)

#------- ice and land:
fname = "landice."+res+"."+tag+".f"+hr+".grib2"
grblandice = pygrib.open(fname)
print("landice = ",grblandice)

for x in grblandice:
  print(x.shortName, x.name, x.level, x.typeOfLevel)
  y = x.values
  if (x.shortName == "lsm"):
    land = y
  elif (x.shortName == "ci"):
    ice  = y
  print(y.max(), y.min())

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
print("npts = ",npts, flush=True)

#----------------------------------------------------
# prepare variables
histogram = np.zeros(8000)
rate_axis = np.zeros(len(histogram))

sum = 0.0
sumsq = 0.0
sumarea = 0.0
icing_plus = 0.0

#frun = open('running_input','rb')
#binary = frun.read()

dtime = int(3)
#-------------------------------------------------------------------
for hour in range(0, 80*dtime + 1, dtime):
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
        #mapper.locate(j, i, ll)
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
  cmap = plt.get_cmap('PiYG')
  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  fig, ax0 = plt.subplots(nrows=1)

  ax0.set(xlabel = "Longitude")
  ax0.set(ylabel = "Latitude")
  ax0.set_title("Icing Rate (cm/hr) -- " + "{:03d}".format(hour) + " hour lead forecast")

  y1 = int(int(ny*3)/4)
  y2 = int(ny-(int(10/abs(dy))))
  im = ax0.pcolormesh(lons, lats[y1:y2], icing_rate[y1:y2,:], cmap=cmap, norm=norm, alpha = .990)
  fig.colorbar(im, ax=ax0)

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