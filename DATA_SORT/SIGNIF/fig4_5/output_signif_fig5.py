import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import sys

from smyleutils import averaging_utils as avg
from smyleutils import bootstrap_utils as boot

initmon=['11','02','09']

basepath="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/Uzm/"
savdir="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SIGNIF/fig4_5/"

for init in initmon:
    high = xr.open_dataset(basepath+'Uzm_BSMYLE-CW3E-L83_day_init'+init+'.nc')
    low = xr.open_dataset(basepath+'Uzm_BSMYLE-CW3E_day_init'+init+'.nc')
    era5 = xr.open_dataset(basepath+'Uzm_ERA5_day_init'+init+'.nc')

    # select the first 6 months since that's what we have for the high top
    startdate = high.time.isel(time=0).values ; enddate = high.time.isel(time=high.time.size-1).values
    high = high.sel(time=slice(startdate,enddate)).Uzm
    low = low.sel(time=slice(startdate,enddate)).Uzm
    era5 = era5.sel(time=slice(startdate,enddate)).Uzm

    high_10 = high.interp(ilev=10)
    low_10 = low.interp(ilev=10)
    era5_10 = era5.sel(level=10.)

    high_10_nh = avg.cosweightlat(high_10, 60, 70)
    high_10_sh = avg.cosweightlat(high_10, -70, -60)
    low_10_nh = avg.cosweightlat(low_10, 60, 70)
    low_10_sh = avg.cosweightlat(low_10, -70, -60)
    era5_10_nh = avg.cosweightlat(era5_10, 60, 70)
    era5_10_sh = avg.cosweightlat(era5_10, -70, -60)

    high_10_nh_clim = high_10_nh.mean('init_year') ; high_10_sh_clim = high_10_sh.mean('init_year')
    low_10_nh_clim = low_10_nh.mean('init_year') ; low_10_sh_clim = low_10_sh.mean('init_year')
    era5_10_nh_clim = era5_10_nh.mean('init_year') ; era5_10_sh_clim = era5_10_sh.mean('init_year')

    high_10_nh_anoms = high_10_nh - high_10_nh_clim ; high_10_sh_anoms = high_10_sh - high_10_sh_clim
    low_10_nh_anoms = low_10_nh - low_10_nh_clim ; low_10_sh_anoms = low_10_sh - low_10_sh_clim
    era5_10_nh_anoms = era5_10_nh - era5_10_nh_clim ; era5_10_sh_anoms = era5_10_sh - era5_10_sh_clim

    mons = int(init) + np.arange(0,6,1)
    mons = np.where(mons <= 12, mons, mons-12)


    high_10_sh_mon = [ 
     high_10_sh_anoms.where( high_10_sh_anoms.time.dt.month == imon, drop=True).mean('time') for imon in mons]
    high_10_sh_mon = xr.concat(high_10_sh_mon, dim=mons)
    high_10_sh_mon = high_10_sh_mon.rename(concat_dim='month')

    high_10_nh_mon = [
     high_10_nh_anoms.where( high_10_nh_anoms.time.dt.month == imon, drop=True).mean('time') for imon in mons]
    high_10_nh_mon = xr.concat(high_10_nh_mon, dim=mons)
    high_10_nh_mon = high_10_nh_mon.rename(concat_dim='month')

    low_10_sh_mon = [
     low_10_sh_anoms.where( low_10_sh_anoms.time.dt.month == imon, drop=True).mean('time') for imon in mons]
    low_10_sh_mon = xr.concat(low_10_sh_mon, dim=mons)
    low_10_sh_mon = low_10_sh_mon.rename(concat_dim='month')

    low_10_nh_mon = [
     low_10_nh_anoms.where( low_10_nh_anoms.time.dt.month == imon, drop=True).mean('time') for imon in mons]
    low_10_nh_mon = xr.concat(low_10_nh_mon, dim=mons)
    low_10_nh_mon = low_10_nh_mon.rename(concat_dim='month')

    era5_10_sh_mon = [
     era5_10_sh_anoms.where( era5_10_sh_anoms.time.dt.month == imon, drop=True).mean('time') for imon in mons]
    era5_10_sh_mon = xr.concat(era5_10_sh_mon, dim=mons)
    era5_10_sh_mon = era5_10_sh_mon.rename(concat_dim='month')

    era5_10_nh_mon = [
     era5_10_nh_anoms.where( era5_10_nh_anoms.time.dt.month == imon, drop=True).mean('time') for imon in mons]
    era5_10_nh_mon = xr.concat(era5_10_nh_mon, dim=mons)
    era5_10_nh_mon = era5_10_nh_mon.rename(concat_dim='month')

    # first bootstrap with replacement across the member dimension (100 samples)
    high_10_nh_mon = high_10_nh_mon.transpose('M',...)
    high_10_sh_mon = high_10_sh_mon.transpose('M',...)
    low_10_nh_mon = low_10_nh_mon.transpose('M',...)
    low_10_sh_mon = low_10_sh_mon.transpose('M',...)

    boot_high_nh_mem = boot.bootgen(high_10_nh_mon, nboots=100)
    boot_low_nh_mem = boot.bootgen(low_10_nh_mon, nboots=100)

    boot_high_sh_mem = boot.bootgen(high_10_sh_mon, nboots=100)
    boot_low_sh_mem = boot.bootgen(low_10_sh_mon, nboots=100)

    boot_high_nh_mem = boot_high_nh_mem.transpose("init_year",...) # move init year to left most
    boot_high_nh_mem = boot_high_nh_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap dim
    boot_high_nh_mem = boot_high_nh_mem.rename({'isample':'M'}) # renaming the member dim

    boot_low_nh_mem = boot_low_nh_mem.transpose("init_year",...) # move init year to left most
    boot_low_nh_mem = boot_low_nh_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap dim
    boot_low_nh_mem = boot_low_nh_mem.rename({'isample':'M'}) # renaming the member dim

    boot_high_sh_mem = boot_high_sh_mem.transpose("init_year",...) # move init year to left most
    boot_high_sh_mem = boot_high_sh_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap dim
    boot_high_sh_mem = boot_high_sh_mem.rename({'isample':'M'}) # renaming the member dim

    boot_low_sh_mem = boot_low_sh_mem.transpose("init_year",...) # move init year to left most
    boot_low_sh_mem = boot_low_sh_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap dim
    boot_low_sh_mem = boot_low_sh_mem.rename({'isample':'M'}) # renaming the member dim

    # loop over the bootstrap samples from the first round of bootstrapping and 
    # bootstrap over time chunks
    allboots_high_nh_cor=[] ; allboots_high_sh_cor=[]
    allboots_low_nh_cor=[] ; allboots_low_sh_cor=[]
    for iboot in np.arange(0,boot_high_nh_mem.ibootmem.size,1):
        print(iboot)
        boot_high_nh_time = boot.bootgenchunk_multimem(boot_high_nh_mem.isel(ibootmem=iboot),
                             5,11,nboots=24, seed=iboot+1)
        boot_high_nh_time_stack = boot_high_nh_time.stack(init_year=['imem','isample'])
        boot_high_nh_time_stack = boot_high_nh_time_stack.isel(init_year=slice(0,high_10_nh_mon.init_year.size))
        boot_high_nh_time_stack = boot_high_nh_time_stack.mean('M')

        boot_low_nh_time = boot.bootgenchunk_multimem(boot_low_nh_mem.isel(ibootmem=iboot),
                             5,11,nboots=24, seed=iboot+1)
        boot_low_nh_time_stack = boot_low_nh_time.stack(init_year=['imem','isample'])
        boot_low_nh_time_stack = boot_low_nh_time_stack.isel(init_year=slice(0,low_10_nh_mon.init_year.size))
        boot_low_nh_time_stack = boot_low_nh_time_stack.mean('M')

        era5_10_nh_mon = era5_10_nh_mon.transpose("init_year",...)
        boot_obs_nh = boot.bootgenchunk_multimem(era5_10_nh_mon, 5, 11, nboots=100, seed=iboot+1)
        boot_obs_nh_stack = boot_obs_nh.stack(init_year=['imem','isample'])
        boot_obs_nh_stack = boot_obs_nh_stack.isel(init_year=slice(0,era5_10_nh_mon.init_year.size))

        cor_high_nh = xr.corr(boot_obs_nh_stack, boot_high_nh_time_stack, dim='init_year')
        allboots_high_nh_cor.append(cor_high_nh)

        cor_low_nh = xr.corr(boot_obs_nh_stack, boot_low_nh_time_stack, dim='init_year')
        allboots_low_nh_cor.append(cor_low_nh)


        boot_high_sh_time = boot.bootgenchunk_multimem(boot_high_sh_mem.isel(ibootmem=iboot),
                            5,11,nboots=25, seed=iboot+1)
        boot_high_sh_time_stack = boot_high_sh_time.stack(init_year=['imem','isample'])
        boot_high_sh_time_stack = boot_high_sh_time_stack.isel(init_year=slice(0,high_10_sh_mon.init_year.size))
        boot_high_sh_time_stack = boot_high_sh_time_stack.mean('M')

        boot_low_sh_time = boot.bootgenchunk_multimem(boot_low_sh_mem.isel(ibootmem=iboot),
                             5,11,nboots=25, seed=iboot+1)
        boot_low_sh_time_stack = boot_low_sh_time.stack(init_year=['imem','isample'])
        boot_low_sh_time_stack = boot_low_sh_time_stack.isel(init_year=slice(0,low_10_sh_mon.init_year.size))
        boot_low_sh_time_stack = boot_low_sh_time_stack.mean('M')

        era5_10_sh_mon = era5_10_sh_mon.transpose("init_year",...)
        boot_obs_sh = boot.bootgenchunk_multimem(era5_10_sh_mon, 5, 11, nboots=100, seed=iboot+1)
        boot_obs_sh_stack = boot_obs_sh.stack(init_year=['imem','isample'])
        boot_obs_sh_stack = boot_obs_sh_stack.isel(init_year=slice(0,era5_10_sh_mon.init_year.size))

        cor_high_sh = xr.corr(boot_obs_sh_stack, boot_high_sh_time_stack, dim='init_year')
        allboots_high_sh_cor.append(cor_high_sh)

        cor_low_sh = xr.corr(boot_obs_sh_stack, boot_low_sh_time_stack, dim='init_year')
        allboots_low_sh_cor.append(cor_low_sh)

   
    allboots_high_nh_cor = xr.concat(allboots_high_nh_cor, dim='iboot2')
    allboots_high_nh_cor = allboots_high_nh_cor.stack(boot=['iboot','iboot2'])
    min95_high_nh_cor = allboots_high_nh_cor.quantile(0.025, dim='boot')
    min95_high_nh_cor = min95_high_nh_cor.rename('min95_high_nh_cor')
    max95_high_nh_cor = allboots_high_nh_cor.quantile(0.975, dim='boot')     
    max95_high_nh_cor = max95_high_nh_cor.rename('max95_high_nh_cor')

    allboots_high_sh_cor = xr.concat(allboots_high_sh_cor, dim='iboot2')
    allboots_high_sh_cor = allboots_high_sh_cor.stack(boot=['iboot','iboot2'])
    min95_high_sh_cor = allboots_high_sh_cor.quantile(0.025, dim='boot')
    min95_high_sh_cor = min95_high_sh_cor.rename('min95_high_sh_cor')
    max95_high_sh_cor = allboots_high_sh_cor.quantile(0.975, dim='boot')     
    max95_high_sh_cor = max95_high_sh_cor.rename('max95_high_sh_cor')

    allboots_low_nh_cor = xr.concat(allboots_low_nh_cor, dim='iboot2')
    allboots_low_nh_cor = allboots_low_nh_cor.stack(boot=['iboot','iboot2'])
    min95_low_nh_cor = allboots_low_nh_cor.quantile(0.025, dim='boot')
    min95_low_nh_cor = min95_low_nh_cor.rename('min95_low_nh_cor')
    max95_low_nh_cor = allboots_low_nh_cor.quantile(0.975, dim='boot')     
    max95_low_nh_cor = max95_low_nh_cor.rename('max95_low_nh_cor')

    allboots_low_sh_cor = xr.concat(allboots_low_sh_cor, dim='iboot2')
    allboots_low_sh_cor = allboots_low_sh_cor.stack(boot=['iboot','iboot2'])
    min95_low_sh_cor = allboots_low_sh_cor.quantile(0.025, dim='boot')
    min95_low_sh_cor = min95_low_sh_cor.rename('min95_low_sh_cor')
    max95_low_sh_cor = allboots_low_sh_cor.quantile(0.975, dim='boot')     
    max95_low_sh_cor = max95_low_sh_cor.rename('max95_low_sh_cor')

    datout = xr.merge([min95_high_nh_cor, max95_high_nh_cor,
                       min95_high_sh_cor, max95_high_sh_cor,
                       min95_low_nh_cor, max95_low_nh_cor,
                       min95_low_sh_cor, max95_low_sh_cor], compat='override')


    datout.to_netcdf(savdir+'ci_polarvortex_monthly_init'+init+'.nc')
