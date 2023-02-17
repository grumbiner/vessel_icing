import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import ListedColormap

matplotlib.use('Agg') #for batch mode
#matplotlib.use('Qt5Agg') #for interactive mode

import cartopy.crs as ccrs
import cartopy.feature as cfeature

#def icing_plot(icing_rate, lons, lats, left_lon, right_lon, bottom_lat, top_lat, name_tag, date_tag)
 
def icing_plot():
  levels = (0,0.01,0.7,2.0, 4.0, 7.5, 13.5, len(histogram)/10.)
#  cmap = plt.get_cmap('PiYG')
# hand craft a color map
  none = np.array([1., 1., 1., 1.])
  trace = np.array([0., 0., 0., 1.])
  light = np.array([0.5, 0.5, 1., 1])
  medium = np.array([1., 0., 0., 1])
  heavy  = np.array([0., 1., 0., 1])
  extreme = np.array([0., 0., 1., 1])
  vheavy  = np.array([1., 1., 0., 1])
  vvheavy = np.array([0., 1., 1., 1])
  cmap = ListedColormap([none,trace, light, medium, heavy, extreme, vheavy, vvheavy])

  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#cartopy is not good about crossing 180 or 0
# Barents/Kara: 0-90 E, 50-90 N, central_longitude = 45, figure_tag="Barents"
# Okhotsk      120-180 E, 30-80 N

  #proj = ccrs.LambertConformal(central_longitude = 150., central_latitude = 60., cutoff=25.)
  proj = ccrs.PlateCarree()
  ax  = plt.axes(projection = proj)
  fig = plt.figure(figsize=(12,9))
  ax  = fig.add_subplot(1,1,1, projection = proj)

  #ax.set_extent((0, 90, 50, 90), crs = ccrs.PlateCarree() )
  ax.set_extent((115, 165, 35, 65), crs = ccrs.PlateCarree() )
  #ax.set_extent((130,180, 40, 90), crs = ccrs.PlateCarree() )
  #ax.gridlines(crs=ccrs.PlateCarree(),
  #             #xlocs=[180, 190, 200, 220, 240],
  #             xlocs=[150, 160, 170],
  #             ylocs=[40,50,60,66.6,70,80] )
  ax.gridlines(crs=ccrs.PlateCarree() )

  #ax.add_feature(cfeature.GSHHSFeature(levels=[1,2,3,4], scale="l") )
  ax.add_feature(cfeature.GSHHSFeature(levels=[1,2], scale="l") )
  #ax.coastlines()

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
