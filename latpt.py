from math import *

from const import *
#Robert Grumbine
#1 June 2018

#Class for working with latitude - longitude points.
#Future: use haversine formula for distance


class latpt:

  def __init__(self,lat = 0, lon = 0):
    self.lat = lat
    self.lon = lon

  def show(self):
    print (self.lon,self.lat)
  
  def distance(self, x ):
### These are in the const.py
#    earth_radius = 6371.2 #km
#    rpdg = 3.1416 / 180.
##
    ab = (90.-self.lat)*const.rpdg
    ac = (90.-x.lat)   *const.rpdg
    bc = abs(self.lon - x.lon)*const.rpdg

    arg = cos(ab)*cos(ac)+sin(ab)*sin(ac)*cos(bc)
    if (arg > 1):
      arg = 1
    if (arg < -1):
      arg = -1

    return const.earth_radius * acos(arg)
