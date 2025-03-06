import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import sys

from smyleutils import averaging_utils as avg
from smyleutils import qboplot_utils as qbo
from smyleutils import colorbar_utils as cbars
from smyleutils import bootstrap_utils as boot


#initmon=['11','02','09']
#startdate=['1970-12-01','1970-03-01','1970-10-01']
#enddate=['1971-02-28','1970-05-31','1970-12-31']

initmon=['02','09']
startdate=['1970-03-01','1970-10-01']
enddate=['1970-05-31','1970-12-31']

basepath="/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/Uzm/"
savdir="/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SIGNIF/fig3/"

# Mean squared skill score calculation
def calcmsss(mod,obs,dim='init_year'):
    mse_mod = (1./mod[dim].size)*((mod - obs)**2).sum('init_year')
    mse_obs = (1./mod[dim].size)*(obs**2).sum('init_year')
    msss = 1 - (mse_mod / mse_obs)

    # dealing with the levels where low top doesn't have any data
    msss = msss.where( msss != 1, nan)
    return msss


for i in np.arange(0,len(initmon),1):
    high = xr.open_dataset(basepath+'Uzm_BSMYLE-CW3E-L83_day_init'+initmon[i]+'.nc')
    low = xr.open_dataset(basepath+'Uzm_BSMYLE-CW3E_day_init'+initmon[i]+'.nc')
    low['lat'] = high.lat
    era5 = xr.open_dataset(basepath+'Uzm_ERA5_day_init'+initmon[i]+'.nc')
    era5['lat'] = high.lat

    # Pick out the right season
    high = high.sel(time=slice(startdate[i],enddate[i])).mean('time').Uzm
    low = low.sel(time=slice(startdate[i],enddate[i])).mean('time').Uzm
    era5 = era5.sel(time=slice(startdate[i],enddate[i])).mean('time').Uzm

    # interpolate the model data from the pressure levels of the CAM TEM daignostics to ERA5
    high_interp = high.interp(ilev=era5.level)
    low_interp = low.interp(ilev=era5.level)

    # Calculate the ensemble mean
    high_interpm = high_interp.mean('M')
    low_interpm = low_interp.mean('M')

    # Calculate the climatology
    high_clim = high_interpm.mean('init_year')
    low_clim = low_interpm.mean('init_year')
    era5_clim = era5.mean('init_year')

    # Subtract the climatology
    high_interp = high_interp - high_clim
    low_interp = low_interp - low_clim
    era5 = era5 - era5_clim

    # First bootstrap with replacement across the member dimension
    high_interp = high_interp.transpose('M',...)
    low_interp = low_interp.transpose('M',...)
    boot_high_mem = boot.bootgen(high_interp,nboots=100)
    boot_low_mem = boot.bootgen(low_interp,nboots=100)

    # move init_year to the left most dimension for the next round of bootstrapping
    boot_high_mem = boot_high_mem.transpose('init_year',...)
    boot_high_mem = boot_high_mem.rename({'iboot':'ibootmem'})
    boot_high_mem = boot_high_mem.rename({'isample':'M'})

    boot_low_mem = boot_low_mem.transpose('init_year',...)
    boot_low_mem = boot_low_mem.rename({'iboot':'ibootmem'})
    boot_low_mem = boot_low_mem.rename({'isample':'M'})

    era5 = era5.transpose('init_year',...)

    allboots_high_cor=[] ; allboots_high_msss=[]
    allboots_low_cor=[] ; allboots_low_msss=[]
    for iboot in np.arange(0,boot_high_mem.ibootmem.size,1):
        print(iboot)
        boot_high_time = boot.bootgenchunk_multimem(
         boot_high_mem.isel(ibootmem=iboot),5,11,nboots=25,seed=iboot+1)
        boot_high_time_stack = boot_high_time.stack(init_year=['imem','isample'])
        boot_high_time_stack = boot_high_time_stack.isel(init_year=slice(0,high_interp.init_year.size))
        boot_high_time_stack = boot_high_time_stack.mean('M')

        boot_low_time = boot.bootgenchunk_multimem(
         boot_low_mem.isel(ibootmem=iboot),5,11,nboots=25,seed=iboot+1)
        boot_low_time_stack = boot_low_time.stack(init_year=['imem','isample'])
        boot_low_time_stack = boot_low_time_stack.isel(init_year=slice(0,low_interp.init_year.size))
        boot_low_time_stack = boot_low_time_stack.mean('M')

        boot_obs = boot.bootgenchunk_multimem(era5, 5, 11, nboots=100, seed=iboot+1)
        boot_obs_stack = boot_obs.stack(init_year=['imem','isample'])
        boot_obs_stack = boot_obs_stack.isel(init_year=slice(0,era5.init_year.size))

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


    datout.to_netcdf(savdir+'ci_latpre_initmon'+initmon[i]+'.nc')


















   





