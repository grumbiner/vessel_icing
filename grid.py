import sys
from struct import *
import numpy as np

#Includes only the mapping, not the data
#Robert Grumbine
#1 June 2018

import const
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

  def locate(self, i, j, z):
    z.lat = self.firstlat + j*self.dlat
    z.lon = self.firstlon + i*self.dlon
    while (z.lon > 360.):
      z.lon -= 360.

  def cellarea(self, i, j):
    tlat = self.firstlat + j*self.dlat
    dy = (self.dlat) * 111.1
    dx = (self.dlon) * 111.1 * cos(math.pi *tlat / 180.0)
    return abs(dx*dy)

#  def locate(self, i, j):
#    z = latpt()
#    z.lat = self.firstlat + j*self.dlat
#    z.lon = self.firstlon + i*self.dlon
#    while (z.lon > 360.):
#      z.lon -= 360.
#    return z

#############################################################

class global_5min(llgrid):

  def __init__(self, nx = 12*360, ny = 12*180):
    self.dlat = -1./12.
    self.dlon =  1./12.
    self.firstlat = 90 + self.dlat/2. 
    self.firstlon = self.dlon / 2.
    self.nx = nx
    self.ny = ny

class global_halfdeg(llgrid):

  def __init__(self, nx = 720, ny = 360):
    self.dlat = -0.5
    self.dlon =  0.5
    self.firstlat = 90 + self.dlat/2. 
    self.firstlon = self.dlon / 2.
    self.nx = nx
    self.ny = ny

class global_nthdeg(llgrid):

  def __init__(self, n = 1.0):
    self.dlat = -1./float(n)
    self.dlon =  1./float(n)
    self.firstlat = 90.0 + self.dlat/2.
    self.firstlon = self.dlon / 2.
    self.nx = 360*int(n)
    self.ny = 180*int(n)

class global_nmc_nthdeg(llgrid):

  def __init__(self, n = 1):
    self.dlat = -1./float(n)
    self.dlon =  1./float(n)
    self.firstlat = 90.0
    self.firstlon = self.dlon / 2.
    self.nx = 360*int(n)
    self.ny = 180*int(n)+1

