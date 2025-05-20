import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
from CASutils import colormap_utils as mycolors

import cartopy as cart
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from cartopy.feature import NaturalEarthFeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import cartopy.io.shapereader as shpreader

def contourmap_bothoceans_tropics_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, contourlines=False):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree(central_longitude=200))
    #ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    #ax.cmap.set_over(mymap(len(mymap)-1))
    ax.set_aspect('auto')
    ax.set_extent([-180,180,-30,30], crs = ccrs.PlateCarree(central_longitude=200))

    if (labels):
        ax.set_xticks([60-200, 120-200, 180-200, 240-200, 300-200, 360-200])
        ax.set_xticklabels(['60E','120E', '180W', '120W', '60W','0W'], fontsize=fontsize-3)
        ax.set_yticks([-30,-20,-10,0,10,20,30], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['30S','20S','10S','0','10N','20N','30N'], fontsize=fontsize-3)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()


    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="max",
            transform=ccrs.PlateCarree())

    if (contourlines):
        clevs2 = clevs[ np.abs(clevs) > ci/2 ]
        ax.contour(lon,lat,dat,levels=clevs2, colors='black', transform=ccrs.PlateCarree())



    ax.add_feature(cfeature.COASTLINE)

    return ax


def contourmap_bothoceans_robinson_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", fontsize=15, signifdat=None, stipplesignif=False, contourlines=None, contourlinescale=1):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
#    clevs[np.abs(clevs) < ci/2.] = 0
#    print(clevs)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=240))
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)
   # ax.set_extent([-180,180,0,90], crs = ccrs.PlateCarree())

   # if (labels):
        #ax.set_xticks([-180, -120, -60, 0,60,120, 180], crs = ccrs.PlateCarree())
        #ax.set_xticklabels(['180W','120W','60W','0','60E','120E','180E'], fontsize=fontsize-3)
        #ax.set_yticks([0,30,60,90], crs = ccrs.PlateCarree())
        #ax.set_yticklabels(['0','30N','60N','90N'], fontsize=fontsize-3)
        #ax.xformatter = LongitudeFormatter()
        #ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())

    if ( signifdat is not None ):
        lonsignif = signifdat.lon
        signifdat, lonsignif = add_cyclic_point( signifdat, coord=lonsignif)
        if (stipplesignif):
            density=3
            ax.contourf(lonsignif, lat, signifdat, levels=[0,0.5,1], colors='none',
               hatches=[density*'.',density*'.', density*','],
               transform = ccrs.PlateCarree())
        else:
            ax.contourf(lonsignif, lat, signifdat, levels=[0,0.5,1], colors='lightgray',
               transform = ccrs.PlateCarree())

    ax.set_global()

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci]
        ax.contour(lon, lat, dat, levels=clevlines, colors='black', transform=ccrs.PlateCarree())

    return ax




