import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd
import sys
from smyleutils import calendar_utils as cal

dat = xr.open_dataset("/glade/campaign/cgd/cas/observations/sst/ersstv5.185401-202312.nc")
time = cal.YYYYMM2date(dat.time)
dat['time'] = time

ystart=1970 ; yend = 2020

alldat=[]
for iyear in np.arange(ystart,yend+1,1):
    datuse = dat.sel(time=slice(str(iyear)+"-11-01",str(iyear+1)+"-04-30"))
    datuse['time'] = pd.date_range("1970-11","1971-04", freq='MS') + pd.DateOffset(days=15)
    alldat.append(datuse)

alldat = xr.concat(alldat, dim='init_year')
alldat['init_year'] = np.arange(ystart,yend+1,1)
alldat.to_netcdf("/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SST/"+
                 "SST_ERSSTv5_mon_init11.nc")
