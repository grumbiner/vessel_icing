import sys
from struct import *
import numpy as np

#Includes only the mapping, not the data
#Robert Grumbine
#1 June 2018

from const import *
from ijpt import *
from latpt import *

############### Friends to the class
def ok(x, loc):
  return (loc.i >= 0 and loc.i < x.shape[0] and
          loc.j >= 0 and loc.j < x.shape[1] )

def binin(fin, nx, ny):
  binary = fin.read()
  fmt=fmt=str(nx*ny)+'f'
  tmpx = unpack(fmt, binary[0:0+4*nx*ny])
  tmp  = np.zeros((nx,ny))
  count = 0
  for val2 in tmpx:
    j = int(count / nx)
    i = int(count % nx)
    tmp[i,j] = val2
    count += 1
  return tmp

def to_2d(binary, nx, ny, k):
  nb = 4*nx*ny
  fmt=str(nx*ny)+'f'
  tmpx = unpack(fmt, binary[k*nb:(k+1)*nb])
  tmp  = np.zeros((nx,ny))
  count = int(0)
  for val2 in tmpx:
    j = int(count / nx)
    i = int(count % nx)
    tmp[i,j] = val2
    count += 1
  return tmp

#beta is much, much faster:
def to2d_beta(binary, nx, ny, k):
  nb = 4*nx*ny
  fmt = str(nx*ny)+'f'
  tmpx = unpack(fmt, binary[k*nb:(k+1)*nb])
  tmp  = np.asarray(tmpx)
  tmp.shape = (nx, ny)
  return tmp

# Declaring a polar stereographic grid and a latitude-longitude grid
#############################################################
class psgrid:

  def locate(self, i, j, z):
    z.lat = i
    z.lon = j

class llgrid:
  #dlat = float(-1.)
  #dlon = float(1.)
  #firstlat = 90. + dlat/2.
  #firstlon = dlon / 2. 
  #nx = int(360./abs(dlon))
  #ny = int(180./abs(dlat))

  def __init__(self, dlat, dlon, firstlat, firstlon, nx, ny):
    print("hello frmo llgrid.__init__")
    self.dlat         = dlat
    self.dlon         = dlon
    self.firstlat     = firstlat
    self.firstlon     = firstlon
    self.nx           = nx
    self.ny           = ny
    self.darea_base   = abs(self.dlat*self.dlon)*const.degree_area
    self.dlat_rad     = self.dlat * const.rpdg
    self.firstlat_rad = self.firstlat * const.rpdg

  def locate(self, j, i, z):
    z.lat = self.firstlat + j*self.dlat
    z.lon = self.firstlon + i*self.dlon
    while (z.lon > 360.):
      z.lon -= 360.

  def cellarea(self, j, i):
    #Original:
    #tlat = self.firstlat + j*self.dlat
    #dy = (self.dlat) * 111.1
    #dx = (self.dlon) * 111.1 * cos(tlat*const.rpdg )
    #return abs(dx*dy)
    #Promote/pre-compute grid constants (darea_base) and 
    #    pre-translate degrees to radians in the base class (*lat_rad)
    return (self.darea_base * cos(self.firstlat_rad + j*self.dlat_rad) ) 


#############################################################

class global_5min(llgrid):

  def __init__(self):
    llgrid.__init__(self, -1./12., 1/12., 90 - 1./24., 1./24., 12*360, 12*180) 

class global_halfdeg(llgrid):

  def __init__(self, nx = 720, ny = 360):
    llgrid.__init__(self, -0.5, 0.5, 90 - 0.25, 0.25, 720, 360)

class global_nthdeg(llgrid):

  def __init__(self, n = 1.0):
    fn = float(n)
    llgrid.__init__(self, -1./fn, 1./fn, 90.,  1./2./fn, n*360, n*180) 

class global_nmc_nthdeg(llgrid):

  def __init__(self, n = 1):
    fn = float(n)
    llgrid.__init__(self, -1./fn, 1./fn, 90.,  1./2./fn, n*360, n*180+1) 

