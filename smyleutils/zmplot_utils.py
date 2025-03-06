import matplotlib.pyplot as plt
import numpy as np
from CASutils import colormap_utils as  mycolors
import sys

def plotlatlogpre_to10(fig, data, lat, pre, ci, cmin, cmax, titlestr, x1=0.1, x2=0.9, y1=0.1, y2=0.9, ylabel=True):
    """
    Plot a pressure versus latitude contour plot up to 0.01hPa.
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
    mymap = mycolors.blue2red_cmap(nlevs)

    plt.rcParams['font.size'] = '11'

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])

    ax.contourf(lat,-1.*np.log10(pre), data, levels=clevs, cmap=mymap, extend='max')
    ax.contour(lat,-1.*np.log10(pre), data, levels=clevs[ clevs != 0], colors='black', linewidths=0.5)
    ax.set_ylim(-np.log10(1000.),-np.log10(10))
    ax.set_yticks([-np.log10(1000),-np.log10(300),-np.log10(100),-np.log10(30),-np.log10(10)])
    ax.set_yticklabels(['1000','300','100','30','10'])
    if (ylabel):
        ax.set_ylabel('Pressure (hPa)', labelpad=-4)
    ax.set_title(titlestr, fontsize=16)
    ax.set_xlabel('Latitude $^{\circ}$N')

    return ax

def plotlatlogpre_to1(fig, data, lat, pre, ci, cmin, cmax, titlestr, x1=0.1, x2=0.9, y1=0.1, y2=0.9, ylabel=True, fsize=11, yticklabels=True, cmap='blue2red', signifdat=None):
    """
    Plot a pressure versus latitude contour plot up to 0.01hPa.
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    if (cmap == "precip_nowhite"):
        mymap = mycolors.precip_cmap_nowhite(nlevs)



    plt.rcParams['font.size'] = fsize

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1])

    ax.contourf(lat,-1.*np.log10(pre), data, levels=clevs, cmap=mymap, extend='max')

    if (signifdat is not None):
        #print(signifdat)
        #print(lat)
        ax.contourf(lat, -1.*np.log10(pre), signifdat, levels=[0,0.5,1], colors='lightgray')


    ax.contour(lat,-1.*np.log10(pre), data, levels=clevs[ np.abs(clevs) > ci/2], colors='black', linewidths=0.5)
    ax.set_ylim(-np.log10(1000.),-np.log10(1))
    ax.set_yticks([-np.log10(1000),-np.log10(100),-np.log10(10),-np.log10(1)])
    if (yticklabels):
        ax.set_yticklabels(['1000','100','10','1'])
    else:
        ax.set_yticklabels([' ',' ',' ',' '])
    if (ylabel):
        ax.set_ylabel('Pressure (hPa)', labelpad=-4)
    ax.set_title(titlestr, fontsize=fsize+2)
    ax.set_xlabel('Latitude $^{\circ}$N')

    return ax




