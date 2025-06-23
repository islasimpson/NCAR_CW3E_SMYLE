import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import sys

from smyleutils import averaging_utils as avg
from smyleutils import qboplot_utils as qbo
from smyleutils import colorbar_utils as cbars
from smyleutils import bootstrap_utils as boot

# Read in the data
initmon=['11','02','09']
basepath="/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/Uzm/"
savdir="/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SIGNIF/fig1/"

# Mean squared skill score calculation
def calcmsss(mod,obs,dim='init_year'):
    mse_mod = (1./mod[dim].size)*((mod - obs)**2).sum('init_year')
    mse_obs = (1./mod[dim].size)*(obs**2).sum('init_year')
    msss = 1 - (mse_mod / mse_obs)
    
    # dealing with the levels where low top doesn't have any data
    msss = msss.where( msss != 1, nan)
    return msss

for init in initmon:
    high = xr.open_dataset(basepath+'Uzm_BSMYLE-CW3E-L83_day_init'+init+'.nc')
    low = xr.open_dataset(basepath+'Uzm_BSMYLE-CW3E_day_init'+init+'.nc')
    era5 = xr.open_dataset(basepath+'Uzm_ERA5_day_init'+init+'.nc')

    # select the first 6 months since that's what we have for the high top
    startdate = high.time.isel(time=0).values ; enddate = high.time.isel(time=high.time.size-1).values
    high = high.sel(time=slice(startdate,enddate)).Uzm
    low = low.sel(time=slice(startdate,enddate)).Uzm
    era5 = era5.sel(time=slice(startdate,enddate)).Uzm

    high_tr = avg.cosweightlat(high,-5,5).load()
    low_tr = avg.cosweightlat(low,-5,5).load()
    era5_tr = avg.cosweightlat(era5,-5,5).load()

    # interpolate the model data from the pressure levels of the CAM TEM diagnostics to ERA5
    high_tr_interp = high_tr.interp(ilev=era5.level)
    low_tr_interp = low_tr.interp(ilev=era5.level)

    # Calculate the ensemble mean
    high_tr_interpm = high_tr_interp.mean('M')
    low_tr_interpm = low_tr_interp.mean('M')
    
    # calculate the lead dependent climatology
    era5_tr_clim = era5_tr.mean('init_year')
    high_tr_interpm_clim = high_tr_interpm.mean('init_year')
    low_tr_interpm_clim = low_tr_interpm.mean('init_year')

    # subtract the lead dependent climatology
    era5_tr = era5_tr - era5_tr_clim
    high_tr_interp = high_tr_interp - high_tr_interpm_clim
    low_tr_interp = low_tr_interp - low_tr_interpm_clim

    # First bootstrap with replacement across the member dimension (100 samples)
    mod_high = high_tr_interp ; mod_low = low_tr_interp ; obs = era5_tr
    boot_high_mem = boot.bootgen(mod_high, nboots=100)
    boot_low_mem = boot.bootgen(mod_low, nboots=100)

    boot_high_mem = boot_high_mem.transpose("init_year",...) # move init_year to the left most dimension for the next round of bootstrapping
    boot_high_mem = boot_high_mem.rename({'iboot':'ibootmem'}) # remaining the bootstrap sample dimension
    boot_high_mem = boot_high_mem.rename({'isample':'M'}) # renaming the member dimension
    
    boot_low_mem = boot_low_mem.transpose("init_year",...) # move init_year to the left most dimension for the next round of bootstrapping
    boot_low_mem = boot_low_mem.rename({'iboot':'ibootmem'}) # remaining the bootstrap sample dimension
    boot_low_mem = boot_low_mem.rename({'isample':'M'}) # renaming the member dimension

    # loop over the bootstrap samples from the first round to avoid memory issues
    allboots_high_cor=[]
    allboots_high_msss=[]
    allboots_low_cor=[]
    allboots_low_msss=[]
    for iboot in np.arange(0,boot_high_mem.ibootmem.size,1): 
        print(iboot) 
        # ---bootstrap high top
        
        # bootstrap over the init_year dimension using 5 year blocks
        boot_high_time = boot.bootgenchunk_multimem(boot_high_mem.isel(ibootmem=iboot),5,11,nboots=25,seed=iboot+1)
        boot_high_time_stack = boot_high_time.stack(init_year=['imem','isample'])
        boot_high_time_stack = boot_high_time_stack.isel(init_year=slice(0,mod_high.init_year.size))
        # compute ensemble mean
        boot_high_time_stack = boot_high_time_stack.mean('M')
        
        
        # bootstrap over the init_year dimension using 5 year blocks
        boot_low_time = boot.bootgenchunk_multimem(boot_low_mem.isel(ibootmem=iboot),5,11,nboots=25,seed=iboot+1)
        boot_low_time_stack = boot_low_time.stack(init_year=['imem','isample'])
        boot_low_time_stack = boot_low_time_stack.isel(init_year=slice(0,mod_low.init_year.size))
        # compute ensemble mean
        boot_low_time_stack = boot_low_time_stack.mean('M')
        
        
        # bootstrap the obs with the same seed
        boot_obs = boot.bootgenchunk_multimem(obs, 5, 11, nboots=100, seed=iboot+1)
        boot_obs_stack = boot_obs.stack(init_year=['imem','isample'])
        boot_obs_stack = boot_obs_stack.isel(init_year=slice(0,obs.init_year.size))
        
        cor_high = xr.corr(boot_obs_stack, boot_high_time_stack, dim='init_year')
        allboots_high_cor.append(cor_high)
        
        msss_high = calcmsss(boot_high_time_stack, boot_obs_stack)
        allboots_high_msss.append(msss_high)
        
        cor_low = xr.corr(boot_obs_stack, boot_low_time_stack, dim='init_year')
        allboots_low_cor.append(cor_low)
        
        msss_low = calcmsss(boot_low_time_stack, boot_obs_stack)
        allboots_low_msss.append(msss_low)
        

    allboots_high_cor = xr.concat(allboots_high_cor, dim='iboot2')
    allboots_high_msss = xr.concat(allboots_high_msss, dim='iboot2')
    allboots_high_cor = allboots_high_cor.stack(boot=['iboot','iboot2'])
    allboots_high_msss = allboots_high_msss.stack(boot=['iboot','iboot2'])
    
    allboots_low_cor = xr.concat(allboots_low_cor, dim='iboot2')
    allboots_low_msss = xr.concat(allboots_low_msss, dim='iboot2')
    allboots_low_cor = allboots_low_cor.stack(boot=['iboot','iboot2'])
    allboots_low_msss = allboots_low_msss.stack(boot=['iboot','iboot2'])
    
    allboots_dif_cor = allboots_high_cor - allboots_low_cor
    allboots_dif_msss = allboots_high_msss - allboots_low_msss
    
    min95_high_cor = allboots_high_cor.quantile(0.025, dim='boot') ; min95_high_cor = min95_high_cor.rename('min95_high_cor')
    max95_high_cor = allboots_high_cor.quantile(0.975, dim='boot') ; max95_high_cor = max95_high_cor.rename('max95_high_cor')
    
    min95_low_cor = allboots_low_cor.quantile(0.025, dim='boot') ; min95_low_cor = min95_low_cor.rename('min95_low_cor')
    max95_low_cor = allboots_low_cor.quantile(0.975, dim='boot') ; max95_low_cor = max95_low_cor.rename('max95_low_cor')
    
    min95_high_msss = allboots_high_msss.quantile(0.025, dim='boot') ; min95_high_msss = min95_high_msss.rename('min95_high_msss')
    max95_high_msss = allboots_high_msss.quantile(0.975, dim='boot') ; max95_high_msss = max95_high_msss.rename('max95_high_msss')
    
    min95_low_msss = allboots_low_msss.quantile(0.025, dim='boot') ; min95_low_msss = min95_low_msss.rename('min95_low_msss')
    max95_low_msss = allboots_low_msss.quantile(0.975, dim='boot') ; max95_low_msss = max95_low_msss.rename('max95_low_msss')
    
    min95_dif_cor = allboots_dif_cor.quantile(0.025, dim='boot') ; min95_dif_cor = min95_dif_cor.rename('min95_dif_cor')
    max95_dif_cor = allboots_dif_cor.quantile(0.975, dim='boot') ; max95_dif_cor = max95_dif_cor.rename('max95_dif_cor')
    
    min95_dif_msss = allboots_dif_msss.quantile(0.025, dim='boot') ; min95_dif_msss = min95_dif_msss.rename('min95_dif_msss')
    max95_dif_msss = allboots_dif_msss.quantile(0.975, dim='boot') ; max95_dif_msss = max95_dif_msss.rename('max95_dif_msss')
    
    datout = xr.merge([min95_high_cor, max95_high_cor,
                       min95_low_cor, max95_low_cor,
                       min95_high_msss, max95_high_msss,
                       min95_low_msss, max95_low_msss,
                       min95_dif_cor, max95_dif_cor,
                       min95_dif_msss, max95_dif_msss], compat='override')


    datout.to_netcdf(savdir+'ci_pretime_initmon'+init+'.nc')







 
