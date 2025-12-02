import xarray as xr
import pandas as pd
import numpy as np
from pandas import Timedelta as timedelta
import sys

def fixcesmtime_daily(dat, timebndsvar='time_bnds'):
    """ Fix the CESM timestamp using the average of time_bnds
        Add in a fix of the times for the leap years """

    try:
        timebndavg = np.array(dat.isel(M=0)[timebndsvar],
                  dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
    except:
        timebndavg = np.array(dat[timebndsvar],
                  dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')

    dates = pd.DatetimeIndex(timebndavg)
    lyindices = np.argwhere( (dates.month == 2) & (dates.day == 29) )
    if (len(lyindices) > 0):
        for i in lyindices:
            timebndavg[i] = str(dates.year[i].item())+"-02-28T12:00:00"

    dat['time'] = timebndavg

    return dat
