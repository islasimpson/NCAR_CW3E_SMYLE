import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import sys

from smyleutils import averaging_utils as avg
from smyleutils import qboplot_utils as qbo
from smyleutils import colorbar_utils as cbars
from smyleutils import bootstrap_utils as boot


basepath="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/U_400_100_avg/"
savdir="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SIGNIF/fig10/"

def calcmsss(mod,obs,dim='init_year'):
    mse_mod = (1./mod[dim].size)*((mod - obs)**2).sum('init_year')
    mse_obs = (1./mod[dim].size)*(obs**2).sum('init_year')
    msss = 1 - (mse_mod / mse_obs)

    # dealing with the levels where low top doesn't have any data
    msss = msss.where( msss != 1, nan)
    return msss

era5 = xr.open_dataset(basepath+'U_400_100_avg_ERA5_init11.nc').__xarray_dataarray_variable__
l83 = xr.open_dataset(basepath+'U_400_100_avg_L83_init11.nc').__xarray_dataarray_variable__
l32 = xr.open_dataset(basepath+'U_400_100_avg_L32_init11.nc').__xarray_dataarray_variable__

# making sure there are no inconsistencies at high decimal places in the lons and lats
l83['lon'] = era5.lon ; l83['lat'] = era5.lat
l32['lon'] = era5.lon ; l32['lat'] = era5.lat

# calculate DJF
era5_djf = era5.sel(time=slice("1970-12-01","1971-02-28")).mean('time')
l83_djf = l83.sel(time=slice("1970-12-01","1971-02-28")).mean('time')
l32_djf = l32.sel(time=slice("1970-12-01","1971-02-28")).mean('time')

del(era5)
del(l83)
del(l32)

# calculate the ensemble mean
l83_djf_em = l83_djf.mean('M')
l32_djf_em = l32_djf.mean('M')

# calculate the lead dependent climatology
era5_djf_clim = era5_djf.mean('init_year')
l83_djf_clim = l83_djf_em.mean('init_year')
l32_djf_clim = l32_djf_em.mean('init_year')

# subtract the lead dependent climatology
era5_djf = era5_djf - era5_djf_clim
l83_djf = l83_djf - l83_djf_clim
l32_djf = l32_djf - l32_djf_clim

print('before bootstrapping members')
#---First bootstrap with replacement across the member dimension (100 samples)
l83_djf = l83_djf.transpose('M',...)
l32_djf = l32_djf.transpose('M',...)
boot_l83_mem = boot.bootgen(l83_djf, nboots=100)
boot_l32_mem = boot.bootgen(l32_djf, nboots=100)

# move init_year to the left most dimension for the next round of bootstrapping
boot_l83_mem = boot_l83_mem.transpose("init_year",...)
boot_l83_mem = boot_l83_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap sample dimension
boot_l83_mem = boot_l83_mem.rename({'isample':'M'})

boot_l32_mem = boot_l32_mem.transpose("init_year",...)
boot_l32_mem = boot_l32_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap sample dimension
boot_l32_mem = boot_l32_mem.rename({'isample':'M'})

print('after bootstrapping members')

allboots_l83_cor=[]
allboots_l83_msss=[]
allboots_l32_cor=[]
allboots_l32_msss=[]
for iboot in np.arange(0,boot_l83_mem.ibootmem.size,1):
    print(iboot)
    boot_l83_time = boot.bootgenchunk_multimem(
      boot_l83_mem.isel(ibootmem=iboot), 5, 11, nboots=25, seed=iboot+1)
    boot_l83_time_stack = boot_l83_time.stack(init_year=['imem','isample'])
    boot_l83_time_stack = boot_l83_time_stack.isel(init_year=slice(0,l83_djf.init_year.size))
    boot_l83_time_stack = boot_l83_time_stack.mean('M')

    boot_l32_time = boot.bootgenchunk_multimem(
      boot_l32_mem.isel(ibootmem=iboot), 5, 11, nboots=25, seed=iboot+1)
    boot_l32_time_stack = boot_l32_time.stack(init_year=['imem','isample'])
    boot_l32_time_stack = boot_l32_time_stack.isel(init_year=slice(0,l32_djf.init_year.size))
    boot_l32_time_stack = boot_l32_time_stack.mean('M')

    boot_obs = boot.bootgenchunk_multimem(era5_djf, 5, 11, nboots=100, seed=iboot+1)
    boot_obs_stack = boot_obs.stack(init_year=['imem','isample'])
    boot_obs_stack = boot_obs_stack.isel(init_year=slice(0,era5_djf.init_year.size))

    cor_l83 = xr.corr(boot_obs_stack, boot_l83_time_stack, dim='init_year')
    allboots_l83_cor.append(cor_l83)

    msss_l83 = calcmsss(boot_l83_time_stack, boot_obs_stack)
    allboots_l83_msss.append(msss_l83)

    cor_l32 = xr.corr(boot_obs_stack, boot_l32_time_stack, dim='init_year')
    allboots_l32_cor.append(cor_l32)

    msss_l32 = calcmsss(boot_l32_time_stack, boot_obs_stack)
    allboots_l32_msss.append(msss_l32)

allboots_l83_cor = xr.concat(allboots_l83_cor, dim='iboot2')
allboots_l83_msss = xr.concat(allboots_l83_msss, dim='iboot2')
allboots_l83_cor = allboots_l83_cor.stack(boot=['iboot','iboot2'])
allboots_l83_msss = allboots_l83_msss.stack(boot=['iboot','iboot2'])

allboots_l32_cor = xr.concat(allboots_l32_cor, dim='iboot2')
allboots_l32_msss = xr.concat(allboots_l32_msss, dim='iboot2')
allboots_l32_cor = allboots_l32_cor.stack(boot=['iboot','iboot2'])
allboots_l32_msss = allboots_l32_msss.stack(boot=['iboot','iboot2'])

allboots_dif_cor = allboots_l83_cor - allboots_l32_cor
allboots_dif_msss = allboots_l83_msss - allboots_l32_msss

min95_l83_cor = allboots_l83_cor.quantile(0.025, dim='boot') ; min95_l83_cor = min95_l83_cor.rename('min95_l83_cor')
max95_l83_cor = allboots_l83_cor.quantile(0.975, dim='boot') ; max95_l83_cor = max95_l83_cor.rename('max95_l83_cor')

min95_l32_cor = allboots_l32_cor.quantile(0.025, dim='boot') ; min95_l32_cor = min95_l32_cor.rename('min95_l32_cor')
max95_l32_cor = allboots_l32_cor.quantile(0.975, dim='boot') ; max95_l32_cor = max95_l32_cor.rename('max95_l32_cor')

min95_l83_msss = allboots_l83_msss.quantile(0.025, dim='boot') ; min95_l83_msss = min95_l83_msss.rename('min95_l83_msss')
max95_l83_msss = allboots_l83_msss.quantile(0.975, dim='boot') ; max95_l83_msss = max95_l83_msss.rename('max95_l83_msss')

min95_l32_msss = allboots_l32_msss.quantile(0.025, dim='boot') ; min95_l32_msss = min95_l32_msss.rename('min95_l32_msss')
max95_l32_msss = allboots_l32_msss.quantile(0.975, dim='boot') ; max95_l32_msss = max95_l32_msss.rename('max95_l32_msss')

min95_dif_cor = allboots_dif_cor.quantile(0.025, dim='boot') ; min95_dif_cor = min95_dif_cor.rename('min95_dif_cor')
max95_dif_cor = allboots_dif_cor.quantile(0.975, dim='boot') ; max95_dif_cor = max95_dif_cor.rename('max95_dif_cor')

min95_dif_msss = allboots_dif_msss.quantile(0.025, dim='boot') ; min95_dif_msss = min95_dif_msss.rename('min95_dif_msss')
max95_dif_msss = allboots_dif_msss.quantile(0.975, dim='boot') ; max95_dif_msss = max95_dif_msss.rename('max95_dif_msss')

datout = xr.merge([min95_l83_cor, max95_l83_cor,
                   min95_l32_cor, max95_l32_cor,
                   min95_l83_msss, max95_l83_msss,
                   min95_l32_msss, max95_l32_msss,
                   min95_dif_cor, max95_dif_cor,
                   min95_dif_msss, max95_dif_msss], compat='override')

datout.to_netcdf(savdir+'ci_initmon11_DJF.nc')



