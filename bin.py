import os
import time

import numpy as np

# test/dummy program for looking at binaries
# Robert Grumbine 
# 11 October 2018

fin = open('sst','rb')

dt = np.dtype('f4')
nx = 360*4
ny = 180*4+1
x = np.zeros((nx,ny),dtype=dt)

x = np.fromfile(fin,dtype=dt,count=nx*ny)
x.shape = (ny, nx)
x -= 273.15

start=time.time()
for j in range (0,ny):
  for i in range (0,nx):
    print ("alpha ",i,j, x[j,i])
end=time.time()
print ("j,i ",end - start)

start = time.time()
for i in range (0,nx):
  for j in range (0,ny):
    print ("beta ",i,j, x[j,i])
end = time.time()
print ("i,j ",end - start)


