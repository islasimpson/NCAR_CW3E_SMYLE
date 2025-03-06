import xarray as xr
import numpy as np

dat = xr.open_dataset("/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/Uzm/"+
                      "Uzm_ERA5_day_init02.nc")
dat = dat.THzm
dat = dat*(dat.level/1000.)**(2./7.)
dat = dat.rename('Tzm')
dat.to_netcdf("/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/Tzm/"+
                  "Tzm_ERA5_day_init02.nc")
