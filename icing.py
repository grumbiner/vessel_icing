import sys

import math
import numpy as np
from struct import *

from const import *
from ijpt  import *
from latpt import *
from grid  import *

#Vessel Icing constants, cm/hr:
icing_a = 2.73E-02
icing_b = 2.91E-04
icing_c = 1.84E-06

#working variables:
loc = ijpt()

#mapper  = global_nmc_nthdeg(1) 
mapper  = global_nmc_nthdeg(4) 
nx = int(mapper.nx)
ny = int(mapper.ny)

dt = np.dtype('f4')
sst  = np.zeros((nx, ny))
land = np.zeros((nx, ny))
ice  = np.zeros((nx, ny))
u10  = np.zeros((nx, ny))
v10  = np.zeros((nx, ny))
t2m  = np.zeros((nx, ny))
speed  = np.zeros((nx, ny))
icing_rate  = np.zeros((nx, ny))
icing_plus = 0.0
histogram = np.zeros(800)

#read static sst, land, ice
#general use:
fmt=str(nx*ny)+'f'

fin = open('sst','rb')
binary = fin.read()
sst = to_2d(binary, nx, ny, 0)
fin.close()


flice = open('landice','rb')
binary = flice.read()
land = to_2d(binary, nx, ny, 0)
ice  = to_2d(binary, nx, ny, 1)
flice.close()

#convert to C if needed
if (sst.max() > 150) :
   sst -= 273.15

print ("sst, land, ice max ",sst.max(), land.max(), ice.max() )

ll = latpt()

sum = 0.0
sumsq = 0.0
sumarea = 0.0
nb=4*nx*ny

frun = open('running_input','rb')
binary = frun.read()

for hour in range(0, 241, 3):
  print("hour = ",hour)
  tau = int(hour / 3)
  t2m = to_2d(binary, nx, ny, 3*tau+0)
  u10 = to_2d(binary, nx, ny, 3*tau+1)
  v10 = to_2d(binary, nx, ny, 3*tau+2)

# transform units, derive quantities:
  if (t2m.max() > 70):
    t2m -= 273.15
  for j in range (0,ny):
    for i in range (0,nx):
      speed[i,j] = sqrt(u10[i,j]*u10[i,j]+v10[i,j]*v10[i,j])

  tf = -1.7
  for j in range (0,ny):
    for i in range (0,nx):
      if (land[i,j] > 0.5 or ice[i,j] > 0.5 or speed[i,j] > 50. or t2m[i,j] > 0 \
       or t2m[i,j] < -40.0 or sst[i,j] < -1.7 or sst[i,j] > 12.0) :
        icing_rate[i,j] = 0.0
        icing_plus = 0.0
      else :
        PR = speed[i,j]*(-1.7 - t2m[i,j])/(1.+.4*(sst[i,j] + 1.7))
        icing_rate[i,j] = PR*(icing_a + icing_b*PR + icing_c*PR*PR)
        PRplus = (2.+speed[i,j])*(-1.7 - (t2m[i,j]-2.) )/(1.+.4*( max(tf, (sst[i,j]-0.5)) + 1.7))
        icing_plus = PRplus*(icing_a + icing_b*PRplus + icing_c*PRplus*PRplus)
        #dI_dP = icing_a + 2.*icing_b*PR + 3.*icing_c*PR*PR
        #dPdw = (1./speed[i,j])
        #dPdTa = -1./(tf-t2m[i,j])
        #dPdTo =  -0.4/(1+0.4*(sst[i,j] - tf)) 
# cm/hr
      #if (icing_rate[i,j] > 0 ) :
      if (icing_plus > 0 ) :
        mapper.locate(i, j, ll)
        print ("grid ", hour, ll.lat, ll.lon, icing_rate[i,j], icing_plus)
        irate = round(icing_rate[i,j]*10.0)
        if (irate >= 0):
          sum += icing_rate[i,j] * mapper.cellarea(i,j)
          sumsq += icing_rate[i,j]*icing_rate[i,j] * mapper.cellarea(i,j)
          sumarea += mapper.cellarea(i,j)
          histogram[int(irate)] += mapper.cellarea(i,j)
  print ("end of loops for hour = ",hour)
  
     

print ("average, rms nonzero = ",sum / sumarea, sqrt(sumsq/sumarea))
histogram /= sumarea
light    = 0.0 # to 0.7 cm/hr
moderate = 0.0 # to 2.0 cm/hr
heavy    = 0.0 # to 4.0 cm/hr
vheavy   = 0.0
vvheavy  = 0.0
extreme  = 0.0 # over 4.0 cm/hr (std.)

for i in range (0,len(histogram) ):
  rate = float(i) / 10.0
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

print ("low-heavy ",light, moderate, heavy,light+moderate+heavy)
print ("vh-extreme ", vheavy, vvheavy, extreme)
print ("all ",light + moderate + heavy + vheavy + vvheavy + extreme)

#convert_to_c(t2m)
#cm/hr to in/hr, m/s 
#sys.exit()
