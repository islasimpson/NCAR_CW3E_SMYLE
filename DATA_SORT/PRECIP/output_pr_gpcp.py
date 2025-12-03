import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd

dat = xr.open_dataset("/glade/campaign/cgd/cas/observations/pr/gpcp.mon.mean.197901-202412.nc")
time = dat.time.values

ystart=1979 ; yend=2020
alldat=[]
for iyear in np.arange(ystart,yend+1,1):
    datuse = dat.sel(time=slice(str(iyear)+"-11-01",str(iyear+1)+"-04-30"))
    datuse['time'] =pd.date_range("1970-11","1971-04", freq='MS') + pd.DateOffset(days=15)
    alldat.append(datuse)
alldat = xr.concat(alldat, dim='init_year')
alldat['init_year'] = np.arange(ystart,yend+1,1)
alldat.to_netcdf("/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/PRECIP/"+
                      "PRECIP_GPCP_mon_init11.nc")
