import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit
import sys
import pandas as pd

pathout="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/MJO_OLR/"
basepath="/glade/campaign/cgd/cas/observations/JRA3Q/day/OLR/"

dat = xr.open_mfdataset(basepath+"/*.nc").__xarray_dataarray_variable__
dat = dat.where( ~( (dat.time.dt.month == 2) & (dat.time.dt.day == 29)), drop=True)
dat.time.encoding['calendar'] = 'noleap'
dat = dat.sel(time=slice("1979-01-01","2023-12-31"))

ystart=1979
yend=2023
mjofilt=[]
for iyear in np.arange(ystart,yend,1):
    print(iyear)
    datuse = dat.sel(time=slice(str(iyear)+"-11-01",str(iyear+1)+"-04-30"))
    datclim = datuse.mean('time')
    datanoms = datuse - datclim
    timenew = pd.date_range("1970-11-01","1971-04-30")
    datanoms['time'] = timenew
 
    mjofilt.append(filt.wkfilter(datanoms,0.15,1,5,20,100, spd=1))

mjofilt = xr.concat(mjofilt, dim='init_year')
mjofilt['init_year'] = np.arange(ystart,yend,1)
mjofilt = mjofilt.sel(time=slice("1970-12-01","1971-02-28"))
mjofilt = mjofilt.rename('MJO_OLR')

mjofilt.to_netcdf(pathout+"JRA3QfilteredOLR_ERA5_init11.nc")
