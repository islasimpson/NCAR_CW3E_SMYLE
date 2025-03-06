import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd

dat = xr.open_dataset("/project/mojave/observations/OBS-SST/ersstv5.185401-202412.nc")
time = dat.time.values
dat['time'] = pd.to_datetime(time, format='%Y%m')

ystart=1970 ; yend = 2020

alldat=[]
for iyear in np.arange(ystart,yend+1,1):
    datuse = dat.sel(time=slice(str(iyear)+"-12-01",str(iyear+1)+"-02-28"))
    datuse['time'] = pd.date_range("1970-12","1971-02", freq='MS') + pd.DateOffset(days=15)
    alldat.append(datuse)

alldat = xr.concat(alldat, dim='init_year')
alldat['init_year'] = np.arange(ystart,yend+1,1)
alldat.to_netcdf("/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SST/"+
                 "SST_ERSSTv5_mon_init11.nc")
