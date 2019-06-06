# -*- coding: utf-8 -*-
"""
@author: Maria Wind
"""

# import several packages and functions
import numpy as np
from netCDF4 import Dataset
from netCDF4 import date2num
from datetime import datetime, timedelta
import xarray as xr
import os
import pandas as pd

from ICC_OBS_sdm_functions_general import sdm_absolute, sdm_relative

#%% DEFINE SOME VARIABLES

def get_timeseries(y,x,start_y,end_y,mod_dir,var):
    # open dataset    
    tmp = xr.open_dataset(mod_dir)
    # select grid point
    tmp = tmp[var][:,y,x]
    # cut out time slice
    tmp_slice = tmp.sel(time=slice(str(start_y), str(end_y)))
    
    mod_data = tmp_slice.data
    try:
        mod_ts = pd.DatetimeIndex(tmp_slice['time'].data)
    except(TypeError):
        mod_ts = tmp_slice['time']
        
    return mod_data, mod_ts

def sdm_wrap(var, obs_dir, mod_cal_dir, mod_bcp_dir, save_dir, start_f, end_f):
    '''
    This functions applies bias correction (scaled distribution mapping) to
    every grid point and saves a grid point time-series as netCDF.
    
    *****************
    INPUT:
    
    var:        variable name (string) (pr, tasmax, tasmin or rsds)
    obs_dir:    path to observation dataset
                dataset should contain data for years 1981-2010
    mod_cal:    path to dataset of model calibration period (historical)
                dataset should contain data for years 1981-2010
    mod_bcp:    path to dataset of model bias correcting period
        
    start_f, end_f: start and end year of future (bias correction) period 
                to perform bias correction only on the calibration period, 
                set the future period to the calibration period (1981-2010)
                
    !!!!!
        Observation and Model data has to cover the exact same domain!    
    !!!!!
    
    '''  
#%%
    if (var == 'pr') | (var == 'rsds'):
        sdm_type = 'relative'
        lower_threshold = 0.1
        distr = 'gamma'    
    
    elif (var == 'tasmax') | (var == 'tasmin'):
        sdm_type = 'absolute'
        distr = 'normal'
    
    elif var == 'hurs':
        sdm_type = 'absolute'
        distr = 'weibull'
    
    elif var == 'sfcWind':
        sdm_type = 'relative'
        lower_threshold = 0.
        distr = 'weibull'


#%%   
    # define calibration period
    start_cal = 1981
    end_cal = 2010    
    
    # load observation data and select time slice
    obs = xr.open_dataset(obs_dir)
    obs = obs.sel(time=slice(str(start_cal), str(end_cal)))
    obs_calib_ts = obs['time'].data
    obs_calib_ts = pd.DatetimeIndex(obs_calib_ts)
    
    # get size of domain
    len_y = obs[var].shape[1]
    len_x = obs[var].shape[2]
    
    #%% GET DATA FOR BIASCORRECTION

    # Loop over every gridpoint within domain
    for y in range(len_y):
        print(str(y)+' of '+str(len_y))
        for x in range(len_x):
        #Read in the data and timeseries (ts) for  
        
        #1) observations   
            obs_calib = obs[var][:,y,x].data
            # check if observations are NOT nan for this grid point, otherwise jump
            # to next grid point

            if np.all(np.isnan(obs_calib)):
                # print 'no obs data available for x=',x,'y=',y
                continue 
            
        #2) model past (calibration)
            mod_calib, mod_calib_ts = get_timeseries(y,x,start_cal,end_cal,mod_cal_dir,var)
            
            # check if model data for this grid point is NOT nan, otherwise jump
            # to next grid point
            if np.all(np.isnan(mod_calib)):
                #print('no model data available for x=',x,'y=',y)
                continue  
            
        #3) model future (rcp)
            mod_bcperiod, mod_bcperiod_ts = get_timeseries(y,x,start_f,end_f,mod_bcp_dir,var)
            
            
#%%
            if sdm_type == 'absolute':
                bc_vals_bcperiod = sdm_absolute(obs_calib, obs_calib_ts, mod_calib, mod_calib_ts, mod_bcperiod, mod_bcperiod_ts, start_f, end_f, distr = distr)
            elif sdm_type == 'relative':
                bc_vals_bcperiod = sdm_relative(lower_threshold, obs_calib, obs_calib_ts, mod_calib, mod_calib_ts, mod_bcperiod, mod_bcperiod_ts, start_f, end_f, distr = distr)            


            #%%  SAVE BIAS CORRECTED DATA FOR GRIDPOINT (Y,X) TO NETCDF FILE 
            
            savedir = save_dir+'/BIASCORR_TS/'
            if not os.path.exists(savedir):
                os.makedirs(savedir)
            
            savename = var+'_y'+str(y)+'_x'+str(x)+'.nc'
            
            # create netCDF file
            # open new netCDF file in write mode:  
            try:
                dataset = Dataset(savedir+savename,'w',format='NETCDF4_CLASSIC')
            except(IOError):
                os.remove(savedir+savename)
                dataset = Dataset(savedir+savename,'w',format='NETCDF4_CLASSIC')
            
            # create dimensions
            dataset.createDimension('time',None)
            
            # create variables
            times = dataset.createVariable('time','f8',('time',))     
            varis = dataset.createVariable(var, 'f4', ('time'))
            
            # add attributes
            times.units = 'days since 1950-01-01T00:00:00Z'
            times.calendar = 'gregorian'
            # --> change to calendar used in original model when reordering data
            
            if var == 'pr':
                varis.units = 'mm'
                varis.long_name = 'total daily precipitation'
            elif var == 'tasmax':
                varis.units = 'degree_Celsius'
                varis.long_name = 'daily maximum near-surface air temperature'
                varis.standard_name = 'air_temperature'
            elif var == 'tasmin':
                varis.units = 'degree_Celsius'
                varis.long_name = 'daily minimum near-surface air temperature'
                varis.standard_name = 'air_temperature'
            elif var == 'rsds':
                varis.units = 'W m-2'
                varis.long_name = 'surface downwelling shortwave flux'
                varis.standard_name = 'surface_downwelling_shortwave_flux_in_air' 
            elif var == 'hurs':
                varis.units = 'percent'
                varis.long_name = 'daily mean relative humidity'
                varis.standard_name = 'relative_humidity'       
            elif var == 'sfcWind':
                varis.units = 'm s-1'
                varis.long_name = '10m wind speed'
                varis.standard_name = 'wind_speed'                
                           
            # write data to netCDF variable
            varis[:] = bc_vals_bcperiod
            
            # fill in times
            dates = [datetime(start_f,1,1)+k*timedelta(days=1) for k in range(bc_vals_bcperiod.shape[0])]
            times[:] = date2num(dates, units=times.units, calendar=times.calendar)
            
            # global attributes
            dataset.title = 'Bias Corrected Model Data for Grid Point y='+str(y)+', x='+str(x)
            dataset.method = 'Bias Correction Method: Scaled Distribution Mapping (Switanek et al., 2017, doi.org/10.5194/hess-21-2649-2017)'
            dataset.source = 'Created with the ICC-OBS Tool (Institute of Meteorology, University of Natural Resources and Life Sciences, Vienna, Austria)'
            
            dataset.close()            
        #end for x        
    #end for y
            
    print('Bias Correction finished!')
    return 
