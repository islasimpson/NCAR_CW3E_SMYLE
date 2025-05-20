import matplotlib.pyplot as plt
import matplotlib as mpl
from smyleutils import colormap_utils as mycolors
import importlib
importlib.reload(mycolors)
import numpy as np


def plotcolorbar(fig, ci, cmin, cmax, titlestr, x1, x2, y1, y2,
   cmap='blue2red', orient='horizontal', posneg='both', ticks=None, fsize=14, nowhite=False,
   contourlines=False, contourlinescale=1):
    """plot a color bar
       Input:
           fig = the figure identified
           ci = the contour interval for the color map
           cmin = the minimum extent of the contour range
           cmax = the maximum extent of the contour range
           titlestr = the label for the color bar
           x1 = the location of the left edge of the color bar
           x2 = the location of the right edge of the color bar
           y1 = the location of the bottom edge of the color bar
           y2 = the location of the top edge of the color bar
           cmap = the color map to be used (only set up for blue2red at the moment)
           orient = the orientation (horizontal or vertical)
           posneg = if "both", both positive and negative sides are plotted
                    if "pos", only the positive side is plotted
                    if "net", only the negative side is plotted
           ticks = user specified ticklabels
           fsize = user specified font size
           contourlines = used to overplot contour lines
           contourlinescale = scale factor for contour lines to be overplotted
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    #clevs = ci * np.arange(cmin/ci, (cmax+ci)/ci, 1) 
    clevs = cmin + np.arange(0,nlevs-1,1)*ci + ci/2.


    #!!! is this right
    #clevs = clevs-ci/2.

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs, nowhite)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs, nowhite)

    if (cmap == "precip_nowhite"):
        mymap = mycolors.precip_cmap_nowhite(nlevs)


    if (cmap == 'red2blue'):
        mymap = mycolors.red2blue_cmap(nlevs, nowhite)

    clevplot=clevs
    if (posneg == "pos"):
        clevplot = clevs[clevs >= 0]
    if (posneg == "neg"):
        clevplot = clevs[clevs <= 0]

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)


    if (ticks):
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
           orientation=orient, norm=norm, values=clevplot, ticks=ticks)
    else:
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
           orientation=orient, norm=norm, values=clevplot)

    clb.ax.tick_params(labelsize=fsize)
    clb.set_label(titlestr, fontsize=fsize+2)

    if (contourlines):
        #clevlines = (clevs-ci/2.)*contourlinescale
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci/2.]
        if (orient=='horizontal'):
            ax.vlines(clevlines[clevlines > 0]-ci/2.,-5,5, colors='black', linestyle='solid')
            ax.vlines(clevlines[clevlines < 0]+ci/2.,-5,5, colors='black', linestyle='dashed')
        if (orient=='vertical'):
            ax.hlines(clevlines[clevlines > 0]-ci/2.,-10,15, colors='black', linestyle='solid')
            ax.hlines(clevlines[clevlines < 0]+ci/2.,-10,15, colors='black', linestyle='dashed')


    return ax


def plotcolorbar_sayc(fig, clevs, titlestr, x1, x2, y1, y2,
   cmap='blue2red', orient='horizontal', posneg='both', ticks=None, ticklabels=None, fsize=14, nowhite=False,
   contourlines=False, contourlinescale=1, nwhite=2):
    """plot a color bar
       Input:
           fig = the figure identified
           titlestr = the label for the color bar
           clevs = the contour levels
           x1 = the location of the left edge of the color bar
           x2 = the location of the right edge of the color bar
           y1 = the location of the bottom edge of the color bar
           y2 = the location of the top edge of the color bar
           cmap = the color map to be used (only set up for blue2red at the moment)
           orient = the orientation (horizontal or vertical)
           posneg = if "both", both positive and negative sides are plotted
                    if "pos", only the positive side is plotted
                    if "net", only the negative side is plotted
           ticks = user specified ticklabels
           fsize = user specified font size
           contourlines = used to overplot contour lines
           contourlinescale = scale factor for contour lines to be overplotted
           nwhite = the number of white contours to have at the beginning of the 
                    Wheeler and Kiladis color bar
    """

    nlevs = len(clevs)
    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs, nowhite)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs, nowhite)

    if (cmap == 'red2blue'):
        mymap = mycolors.red2blue_cmap(nlevs, nowhite)

    if (cmap == 'wk'):
        mymap = mycolors.wkcmap(nlevs, nwhite)

    clevplot=clevs
    if (posneg == "pos"):
        clevplot = clevs[clevs >= 0]
    if (posneg == "neg"):
        clevplot = clevs[clevs <= 0]

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])

#    ci = (np.max(clevs) - np.min(clevs))/len(clevs)
    cvals = np.linspace(-1,1,len(clevs))
    ci = cvals[1] - cvals[0]

#    norm = mpl.colors.Normalize(vmin=clevs[0]+ci/2., vmax = clevs[0] + (len(clevs)-1)*ci -ci/2.)
    norm = mpl.colors.Normalize(vmin=-1 + ci/2., vmax=-1 + (len(clevs)-1)*ci - ci/2.)

#    clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
#        orientation=orient, norm=norm, values=clevs[0] + np.arange(0,len(clevs)-1,1)*ci+ ci/2.,
#        ticks=clevs[0] + np.arange(0,len(clevs),1)*ci)

    if ticks is not None:
        values2 = -1 + np.arange(0,len(clevs),1)*ci
        ticksplot = [ values2[np.array(clevs) == i][0] for i in ticks ]
        ticksplot = list(ticksplot)
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
            orientation=orient, norm=norm, values=-1 + np.arange(0,len(clevs)-1,1)*ci+ ci/2.,
            ticks=ticksplot)
        if ticklabels is not None:
            clb.set_ticklabels(ticklabels)
        else:
            clb.set_ticklabels(ticks)
    else:
        clb = mpl.colorbar.ColorbarBase(ax, cmap=mymap,
            orientation=orient, norm=norm, values=-1 + np.arange(0,len(clevs)-1,1)*ci+ ci/2.,
            ticks=-1 + np.arange(0,len(clevs),1)*ci)
        if ticklabels is not None:
            clb.set_ticklabels(ticklabels)
        else:
            clb.set_ticklabels(clevs)

    clb.ax.tick_params(labelsize=fsize)
    clb.set_label(titlestr, fontsize=fsize+2)

    if (contourlines):
        #clevlines = (clevs-ci/2.)*contourlinescale
        #clevlines = clevs*contourlinescale
        #clevlines = clevlines[np.abs(clevlines) > ci/2.]
        clevlines = -1+np.arange(0,len(clevs),1)*ci
        clevlines = clevlines[ np.abs(clevlines) > ci/2. ]
        if (orient=='horizontal'):
            ax.vlines(clevlines[clevlines >= 0],-5,5, colors='black', linestyle='solid')
            ax.vlines(clevlines[clevlines < 0],-5,5, colors='black', linestyle='dashed')
        if (orient=='vertical'):
            ax.hlines(clevlines[clevlines >= 0],-10,15, colors='black', linestyle='solid')
            ax.hlines(clevlines[clevlines < 0],-10,15, colors='black', linestyle='dashed')


    return ax

