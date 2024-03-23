import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan

from smyleutils import colormap_utils as mycolors

def plot_lev_time(fig, data, time, pre, ci, cmin, cmax, titlestr, x1=None, x2=None, y1=None, y2=None,
    xlim=None, ylim=None, plevvar='ilev', contourlines=False, contourlinescale=1, xlabel='Time (days)'):
    """ 
    Plot a pressure versus latitude time series in log-pressure coordinates
    """

    data = data.transpose(plevvar,"time")

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)
    mymap = mycolors.blue2red_cmap(nlevs)

    plt.rcParams['font.size'] = '12'

    if (x1):
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1])
    else:
        ax = fig.add_axes()

    ax.contourf(time, -1.*np.log10(pre), data, levels=clevs, cmap=mymap, extend='both')
    ax.set_ylim(-np.log10(1000.),-np.log10(1))
    ax.set_yticks([-np.log10(1000),-np.log10(300),-np.log10(100),-np.log10(30),
                  -np.log10(10),-np.log10(3),-np.log10(1)])
    ax.set_yticklabels(['1000','300','100','30','10','3','1'])
    ax.set_ylabel('Pressure (hPa)')
    ax.set_title(titlestr, fontsize=16)
    ax.set_xlabel(xlabel)

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci ]
        ax.contour(time, -1.*np.log10(pre), data, levels=clevlines, colors='black')

    return ax




