import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import sys
from smyleutils import bootstrap_utils as boot

initmon=['09','11','02']

basepath="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/NAO/"
savdir="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SIGNIF/fig8/"

def calcmsss(mod,obs,dim='init_year'):
    mse_mod = (1./mod[dim].size)*((mod - obs)**2).sum('init_year')
    mse_obs = (1./mod[dim].size)*(obs**2).sum('init_year')
    msss = 1 - (mse_mod / mse_obs)
    
    # dealing with the levels where low top doesn't have any data
    msss = msss.where( msss != 1, nan)
    return msss

for init in initmon:
    high = xr.open_dataset(basepath+'NAO_stationbased_L83_initmon'+init+'.nc')
    low = xr.open_dataset(basepath+'NAO_stationbased_L32_initmon'+init+'.nc')
    era5 = xr.open_dataset(basepath+'NAO_stationbased_ERA5_initmon'+init+'.nc')
    
    # select the first 6 months since that's what we have for the high top
    startdate = high.time.isel(time=0).values ; enddate = high.time.isel(time=high.time.size-1).values
    high = high.sel(time=slice(startdate,enddate)).NAO
    low = low.sel(time=slice(startdate,enddate)).NAO
    era5 = era5.sel(time=slice(startdate,enddate)).NAO

    # Calculate 3 momnth rolling seasonal averages
    high_rolling = high.rolling(time=3, center=True, min_periods=3).mean().dropna('time')
    low_rolling = low.rolling(time=3, center=True, min_periods=3).mean().dropna('time')
    era5_rolling = era5.rolling(time=3, center=True, min_periods=3).mean().dropna('time')

    # remove the lead dependent climatology
    high_rolling = high_rolling - high_rolling.mean('init_year')
    low_rolling = low_rolling - low_rolling.mean('init_year')
    era5_rolling = era5_rolling - era5_rolling.mean('init_year')

    # First bootstrap with replacement across the member dimension
    high_rolling = high_rolling.transpose('M',...)
    low_rolling = low_rolling.transpose('M',...)

    high_mem = boot.bootgen(high_rolling, nboots=100)
    low_mem = boot.bootgen(low_rolling, nboots=100)

    high_mem = high_mem.transpose("init_year",...) # move the init year to the left most
    high_mem = high_mem.rename({'iboot':'ibootmem'}) # rename the bootstrap dimension
    high_mem = high_mem.rename({'isample':'M'}) # rename the member dimension

    low_mem = low_mem.transpose("init_year",...) # move the init year to the left most
    low_mem = low_mem.rename({'iboot':'ibootmem'}) # rename the bootstrap dimension
    low_mem = low_mem.rename({'isample':'M'}) # rename the member dimension

    era5_rolling = era5_rolling.transpose("init_year",...)

    # Loop over the bootstrap samples from the first round and botostrap over time chunks
    allboots_high_cor=[] ; allboots_low_cor=[]
    allboots_high_msss=[] ; allboots_low_msss=[]
    for iboot in np.arange(0,high_mem.ibootmem.size,1):
        print(iboot)
        boot_high_time = boot.bootgenchunk_multimem(high_mem.isel(ibootmem=iboot),
                                  5, 11, nboots=25, seed=iboot+1)
        boot_high_time = boot_high_time.stack(init_year=['imem','isample'])
        boot_high_time = boot_high_time.isel(init_year=slice(0,high.init_year.size))
        boot_high_time = boot_high_time.mean('M')
 
        boot_low_time = boot.bootgenchunk_multimem(low_mem.isel(ibootmem=iboot),
                                  5, 11, nboots=25, seed=iboot+1)
        boot_low_time = boot_low_time.stack(init_year=['imem','isample'])
        boot_low_time = boot_low_time.isel(init_year=slice(0,low.init_year.size))
        boot_low_time = boot_low_time.mean('M')

        boot_era5_time = boot.bootgenchunk_multimem(era5_rolling, 5, 11, nboots=25, seed=iboot+1)
        boot_era5_time = boot_era5_time.stack(init_year=['imem','isample'])
        boot_era5_time = boot_era5_time.isel(init_year=slice(0,era5.init_year.size))

        cor_high = xr.corr(boot_high_time, boot_era5_time, dim='init_year')
        allboots_high_cor.append(cor_high)

        cor_low = xr.corr(boot_low_time, boot_era5_time, dim='init_year')
        allboots_low_cor.append(cor_low)

        msss_high = calcmsss(boot_high_time, boot_era5_time)
        allboots_high_msss.append(msss_high)

        msss_low = calcmsss(boot_low_time, boot_era5_time)
        allboots_low_msss.append(msss_low)

    allboots_high_cor = xr.concat(allboots_high_cor, dim='iboot2')
    allboots_high_cor = allboots_high_cor.stack(boot=['iboot','iboot2'])
    min95_high_cor = allboots_high_cor.quantile(0.025, dim='boot')
    min95_high_cor = min95_high_cor.rename('min95_high_cor')
    max95_high_cor = allboots_high_cor.quantile(0.975, dim='boot')
    max95_high_cor = max95_high_cor.rename('max95_high_cor')

    allboots_low_cor = xr.concat(allboots_low_cor, dim='iboot2')
    allboots_low_cor = allboots_low_cor.stack(boot=['iboot','iboot2'])
    min95_low_cor = allboots_low_cor.quantile(0.025, dim='boot')
    min95_low_cor = min95_low_cor.rename('min95_low_cor')
    max95_low_cor = allboots_low_cor.quantile(0.975, dim='boot')
    max95_low_cor = max95_low_cor.rename('max95_low_cor')

    allboots_high_msss = xr.concat(allboots_high_msss, dim='iboot2')
    allboots_high_msss = allboots_high_msss.stack(boot=['iboot','iboot2'])
    min95_high_msss = allboots_high_msss.quantile(0.025, dim='boot')
    min95_high_msss = min95_high_msss.rename('min95_high_msss')
    max95_high_msss = allboots_high_msss.quantile(0.975, dim='boot')
    max95_high_msss = max95_high_msss.rename('max95_high_msss')

    allboots_low_msss = xr.concat(allboots_low_msss, dim='iboot2')
    allboots_low_msss = allboots_low_msss.stack(boot=['iboot','iboot2'])
    min95_low_msss = allboots_low_msss.quantile(0.025, dim='boot')
    min95_low_msss = min95_low_msss.rename('min95_low_msss')
    max95_low_msss = allboots_low_msss.quantile(0.975, dim='boot')
    max95_low_msss = max95_low_msss.rename('max95_low_msss')



    datout = xr.merge([min95_high_cor, max95_high_cor,
                       min95_low_cor, max95_low_cor,
                       min95_high_msss, max95_high_msss,
                       min95_low_msss, max95_low_msss], compat='override')

    datout.to_netcdf(savdir+'ci_NAO_init'+init+'.nc')


