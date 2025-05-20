import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import sys

from smyleutils import averaging_utils as avg
from smyleutils import qboplot_utils as qbo
from smyleutils import colorbar_utils as cbars
from smyleutils import bootstrap_utils as boot
from CASutils import regrid_utils as regrid

basepath="/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/PRECIP/"
savdir="/project/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SIGNIF/fig_precip/"

def calcmsss(mod,obs,dim='init_year'):
    mse_mod = (1./mod[dim].size)*((mod - obs)**2).sum('init_year')
    mse_obs = (1./mod[dim].size)*(obs**2).sum('init_year')
    msss = 1 - (mse_mod / mse_obs)

    # dealing with the levels where low top doesn't have any data
    msss = msss.where( msss != 1, nan)
    return msss

gpcp = xr.open_dataset(basepath+"PRECIP_GPCP_mon_init11.nc").precip

precl = xr.open_dataset(basepath+'PRECL_BSMYLE-CW3E-L83_mon_init11.nc')
precc = xr.open_dataset(basepath+'PRECC_BSMYLE-CW3E-L83_mon_init11.nc')
high = (precl.PRECL + precc.PRECC)*86400.*1000.

precl = xr.open_dataset(basepath+'PRECL_BSMYLE-CW3E-L32_mon_init11.nc')
precc = xr.open_dataset(basepath+'PRECC_BSMYLE-CW3E-L32_mon_init11.nc')
low = (precl.PRECL + precc.PRECC)*86400.*1000.

gpcp_djf = gpcp.sel(time=slice("1970-12-01","1971-02-28")).mean('time')
high_djf = high.sel(time=slice("1970-12-01","1971-02-28")).mean('time')
low_djf = low.sel(time=slice("1970-12-01","1971-02-28")).mean('time')

high_djf = regrid.regrid_conservative(high_djf, high_djf.lon, high_djf.lat,
                 gpcp_djf.lon, gpcp_djf.lat, reuse_wgts=False, 
                 wgtfile="/project/cas/islas/temp/wgtfile_gpcp2.nc")
low_djf = regrid.regrid_conservative(low_djf, low_djf.lon, low_djf.lat,
                 gpcp_djf.lon, gpcp_djf.lat, reuse_wgts=False,
                 wgtfile="/project/cas/islas/temp/wgtfile_gpcp2.nc")
high_djf = high_djf.sel(init_year=slice(1979,2020))
low_djf = low_djf.sel(init_year=slice(1979,2020))


#---calculate the ensemble mean
high_djf_em = high_djf.mean('M')
low_djf_em = low_djf.mean('M')

#---calculate the lead dependent climatology
gpcp_djf_clim = gpcp_djf.mean('init_year')
high_djf_clim = high_djf_em.mean('init_year')
low_djf_clim = low_djf_em.mean('init_year')

#---subtract the lead dependent climatology
gpcp_djf = gpcp_djf - gpcp_djf_clim 
high_djf = high_djf - high_djf_clim
low_djf = low_djf - low_djf_clim

#---First bootstrap with replacement across the member dimension (100 samples)
high_djf = high_djf.transpose("M",...)
low_djf = low_djf.transpose("M",...)
boot_high_mem = boot.bootgen(high_djf, nboots=100)
boot_low_mem = boot.bootgen(low_djf, nboots=100)

#---move init_year to the left most dimension for the next round of bootstrapping
boot_high_mem = boot_high_mem.transpose("init_year",...)
boot_high_mem = boot_high_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap sample dimension
boot_high_mem = boot_high_mem.rename({'isample':'M'})

boot_low_mem = boot_low_mem.transpose("init_year",...)
boot_low_mem = boot_low_mem.rename({'iboot':'ibootmem'}) # renaming the bootstrap sample dimension
boot_low_mem = boot_low_mem.rename({'isample':'M'})


allboots_high_cor=[]
allboots_high_msss=[]
allboots_low_cor=[]
allboots_low_msss=[]
for iboot in np.arange(0,boot_high_mem.ibootmem.size,1):
    print(iboot)
    boot_high_time = boot.bootgenchunk_multimem(
      boot_high_mem.isel(ibootmem=iboot), 5, 11, nboots=25, seed=iboot+1)
    boot_high_time_stack = boot_high_time.stack(init_year=['imem','isample'])
    boot_high_time_stack = boot_high_time_stack.isel(init_year=slice(0,high_djf.init_year.size))
    boot_high_time_stack = boot_high_time_stack.mean('M')

    boot_low_time = boot.bootgenchunk_multimem(
      boot_low_mem.isel(ibootmem=iboot), 5, 11, nboots=25, seed=iboot+1)
    boot_low_time_stack = boot_low_time.stack(init_year=['imem','isample'])
    boot_low_time_stack = boot_low_time_stack.isel(init_year=slice(0,low_djf.init_year.size))
    boot_low_time_stack = boot_low_time_stack.mean('M')

    boot_obs = boot.bootgenchunk_multimem(gpcp_djf, 5, 11, nboots=100, seed=iboot+1)
    boot_obs_stack = boot_obs.stack(init_year=['imem','isample'])
    boot_obs_stack = boot_obs_stack.isel(init_year=slice(0,gpcp_djf.init_year.size))

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


datout.to_netcdf(savdir+'ci_initmon11_DJF.nc')



