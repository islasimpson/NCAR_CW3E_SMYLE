import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan 
from smyleutils import regrid_utils as regrid
from smyleutils import mapplot_utils as mymaps
from smyleutils import averaging_utils as avg
from smyleutils import linfit_utils as linfit
from smyleutils import plothisto_utils as histo
from smyleutils import bootstrap_utils as boot
from smyleutils import colorbar_utils as cbars
import cartopy.crs as ccrs
import warnings
warnings.filterwarnings('ignore')
import sys

savdir="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/SIGNIF/fig14/"

def calcmsss(mod,obs,dim='init_year'):
    mse_mod = (1./mod[dim].size)*((mod - obs)**2).sum('init_year')
    mse_obs = (1./mod[dim].size)*(obs**2).sum('init_year')
    msss = 1 - (mse_mod / mse_obs)

    # dealing with the levels where low top doesn't have any data
    msss = msss.where( msss != 1, nan)
    return msss

# read in qbo for compositing
basepath="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/Uzm/"
uzm_era5_nov=xr.open_dataset(basepath+'Uzm_ERA5_day_init11.nc').Uzm
uzm_era5_nov = uzm_era5_nov.sel(init_year=slice(1979,2020))
uzm_era5_nov_tr = avg.cosweightlat(uzm_era5_nov,-5,5).load()
uqbo = uzm_era5_nov_tr.interp(level=50.).sel(time=slice("1970-12-01","1971-02-28")).mean('time')
uqbo = uqbo - uqbo.mean('init_year')
uqbo_std = uqbo.std('init_year')

# read in precipitation data and calculate the DJF average
basepath="/glade/campaign/cgd/cas/islas/python_savs/NCAR_CW3E_SMYLE/DATA_SORT/PRECIP/"
gpcp = xr.open_dataset(basepath+"PRECIP_GPCP_mon_init11.nc").precip
l83 = xr.open_dataset(basepath+'PRECT_BSMYLE-CW3E-L83_mon_init11.nc').PRECT*86400.*1000.
l32 = xr.open_dataset(basepath+'PRECT_BSMYLE-CW3E-L32_mon_init11.nc').PRECT*86400.*1000.

gpcp = gpcp.sel(time=slice("1970-12-01","1971-02-28")).mean('time')
l83 = l83.sel(time=slice("1970-12-01","1971-02-28")).mean('time')
l32 = l32.sel(time=slice("1970-12-01","1971-02-28")).mean('time')

gpcp = gpcp.sel(init_year=slice(1979,2020))
l83 = l83.sel(init_year=slice(1979,2020))
l32 = l32.sel(init_year=slice(1979,2020))

# remove the lead dependent climatology
gpcp = gpcp - gpcp.mean('init_year')
l83 = l83 - l83.mean('M').mean('init_year')
l32 = l32 - l32.mean('M').mean('init_year')

# regrid the model conservatively onto the GPCP grid
l83 = regrid.regrid_conservative(l83, l83.lon, l83.lat, gpcp.lon, gpcp.lat, reuse_wgts=False,
                                 wgtfile='/glade/derecho/scratch/islas/temp/wgt.nc')
l32 = regrid.regrid_conservative(l32, l32.lon, l32.lat, gpcp.lon, gpcp.lat, reuse_wgts=True,
                                 wgtfile='/glade/derecho/scratch/islas/temp/wgt.nc')


# first bootstrap with replacement across members
l83 = l83.transpose('M',...)
l32 = l32.transpose('M',...)
boot_l83_mem = boot.bootgen(l83, nboots=100)
boot_l32_mem = boot.bootgen(l32, nboots=100)

# move init_year to the left most dimension for the next round of bootstrapping
boot_l83_mem = boot_l83_mem.transpose('init_year',...)
boot_l83_mem = boot_l83_mem.rename({'iboot':'ibootmem'})
boot_l83_mem = boot_l83_mem.rename({'isample':'M'})

boot_l32_mem = boot_l32_mem.transpose('init_year',...)
boot_l32_mem = boot_l32_mem.rename({'iboot':'ibootmem'})
boot_l32_mem = boot_l32_mem.rename({'isample':'M'})

allboots_l83_cor=[]
allboots_l83_msss=[]
allboots_l32_cor=[]
allboots_l32_msss=[]
for iboot in np.arange(0,boot_l83_mem.ibootmem.size,1):
    print(iboot)
    boot_l83_time = boot.bootgenchunk_multimem(
      boot_l83_mem.isel(ibootmem=iboot), 5, 11, nboots=25, seed=iboot+1)
    boot_l83_time_stack = boot_l83_time.stack(init_year=['imem','isample'])
    boot_l83_time_stack = boot_l83_time_stack.isel(init_year=slice(0,l83.init_year.size))
    boot_l83_time_stack = boot_l83_time_stack.mean('M')

    boot_l32_time = boot.bootgenchunk_multimem(
      boot_l32_mem.isel(ibootmem=iboot), 5, 11, nboots=25, seed=iboot+1)
    boot_l32_time_stack = boot_l32_time.stack(init_year=['imem','isample'])
    boot_l32_time_stack = boot_l32_time_stack.isel(init_year=slice(0,l32.init_year.size))
    boot_l32_time_stack = boot_l32_time_stack.mean('M')

    boot_obs = boot.bootgenchunk_multimem(gpcp, 5, 11, nboots=100, seed=iboot+1)
    boot_obs_stack = boot_obs.stack(init_year=['imem','isample'])
    boot_obs_stack = boot_obs_stack.isel(init_year=slice(0,gpcp.init_year.size))

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

datout.to_netcdf(savdir+'pr_ci_initmon11_DJF.nc')
