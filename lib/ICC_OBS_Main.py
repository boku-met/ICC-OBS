"""
Created on Tue Nov 20 16:31:41 2018
@author: mariaw
"""
import glob
import xarray as xr
import os
import datetime
import numpy as np
import sys
import warnings

from ICC_OBS_Functions import cut_domain, write_netcdf, write_netcdf_model, plausibility_check, check_domain, get_ncattrs, plot_result
from ICC_OBS_read_stationdata import read_stationdata
from ICC_OBS_interpolation import merge_no_grad, merge_hgt_grad
from ICC_OBS_bias_correction_general import sdm_wrap
from ICC_OBS_reorder_bcdata import reorder_data
from itertools import combinations


#%%
class dist_angle_conv:
    def __init__(self, value):
        self.earthradius = 6378.
        self.value = value
    
    def to_angle(self):
        self.angle = np.rad2deg(np.arcsin(self.value / self.earthradius))
        return self.angle
    
    def to_dist(self):
        self.dist = self.earthradius * np.sin(self.value*np.pi/180.)
        return self.dist

def distance_great_circle(u,v):
    """ Fast great-circle path distance using Haversine formula
    (at least faster than geopy.distance.distance but less accurate)
    ACHTUNG! Ergebnis muss erst mit Erdradius multipliziert werden !!
    """

    lat_1, lon_1 = u
    lat_2, lon_2 = v
    lat_1, lon_1, lat_2, lon_2 = lat_1*np.pi/180, lon_1*np.pi/180, lat_2*np.pi/180, lon_2*np.pi/180

    d_lon = np.abs(lon_1 - lon_2)

    A = np.power(np.cos(lat_2)*np.sin(d_lon), 2)
    B = np.power(np.cos(lat_1)*np.sin(lat_2) - np.sin(lat_1)*np.cos(lat_2)*np.cos(d_lon), 2)
    C = np.sin(lat_1)*np.sin(lat_2) + np.cos(lat_1)*np.cos(lat_2)*np.cos(d_lon)

    return np.arctan(np.sqrt(A + B) / C)

def main_function(fil_topo, lat_min, lat_max, lon_min, lon_max, param, start_y, end_y, dir_stations, fn_metadata, fil_obsgrid, model_hist, model_rcp, dir_save, dir_save_name,
                  interp_method = None, distance = None, variogram_model = None, neighbor = None):
    '''
    INPUT VARIABLES:

    fil_topo        netcdf file containing the model topography
    lat/lon         latitude and longitude values of bounding box (in decimal degrees)
    param           name of the variable that should be used (pr, tasmax, tasmin, rsds, ws10, rh)
    start_y         start year of the correction period
    end_y           end year of the correction period
    dir_stations    directory that contains the station data
    fn_metadata     file containing the station metadata information
    fil_obsgrid     netcdf file containing gridded observations
    model_hist      directory to historical model files or netcdf file of historical model run
    model_rcp       directory to model scenario (rcp) files or netcdf file of model scenario run

    dir_save        directory that should be used to store data created by the tool
    dir_save_name   name for the directory that is created for the run

    advanced options (predefined values):

    interp_method   kriging (ordinary kriging algorighm) or
                    idw (inverse distance)

    options for inverse distance interpolation (predefined)
    
    neighbor        Minimum number of neighbors needed to perform barnes or 
                    cressman interpolation for a point. Default is 3.
                    See also: https://unidata.github.io/MetPy/latest/api/generated/metpy.interpolate.inverse_distance_to_grid.html
    distance        distance (in km) at which the weights should be 0 (default: 50km)

    options for kriging:

    variogram_model gaussian, spherical, exponential, linear
    

    '''
    warnings.filterwarnings("ignore")

# load topography data from file (fil_topo)
    topo = xr.open_dataset(fil_topo)
# check specified model domain (lat_min, lat_max, lon_min, lon_max)
    err_domain = check_domain(topo, lat_min, lat_max, lon_min, lon_max)
    if err_domain != 0:
        print('The specified latitude and longitude values are outside the model domain!')
        sys.exit('Please choose a different domain!')

    # cut out the domain according to the bounding box
    topo = cut_domain(topo, lat_min, lat_max, lon_min, lon_max)
    lat = topo.lat.data
    lon = topo.lon.data

# Get the directory with the station data and the metadata file.
    # Read station data to array (stations have to be within lat/lon domain)
    print('\nReading stationdata...')
    df_st, df_metadata = read_stationdata(param, dir_stations, fn_metadata, 1981, 2010)
    st_data = df_st.values
    st_lons = df_metadata.lon.values
    st_lats = df_metadata.lat.values
    st_height = df_metadata.height.values

    # check plausibility of station data
    print('\nChecking plausibility of station data...')
    err = plausibility_check(param, st_data, st_lons, st_lats, lat_min, lat_max, lon_min, lon_max)
    if err != 0:
        sys.exit('Please check your station data and restart the tool.')
        
# default interpolation settings:
    if distance != None:
        distance = dist_angle_conv(distance).to_angle() # converts distance input from km to degree
    if not neighbor is None:
        neighbor = int(neighbor)
        
    if interp_method == None:
        if param == 'pr':
            interp_method = 'kriging'
            variogram_model = 'gaussian'

        else:
            interp_method = 'idw'
            points = list(zip(st_lats, st_lons))
            distances = [6378. * distance_great_circle(p1, p2) for p1, p2 in combinations(points, 2)]
            distance = dist_angle_conv(sum(distances) / len(distances)).to_angle() # dient als firstguess fuer IDW
            neighbor=3

# Read gridded observation data (fil_obsgrid)
    print('\nReading gridded observation data...')
    grid = xr.open_dataset(fil_obsgrid)

    # cut out the domain defined at the top
    grid = cut_domain(grid, lat_min, lat_max, lon_min, lon_max)

# Create a new folder under the dir_save directory using the dir_save_name and a time stamp
    now = datetime.datetime.now()
    dir_save_run = dir_save+'/'+dir_save_name+'_'+now.strftime('%Y-%m-%d_%H%M')+'/'
    print('\nCreating directory for this run: '+dir_save_run)
    os.mkdir(dir_save_run)

    dir_tmp = dir_save_run+'TMP/'
    os.mkdir(dir_tmp)
    # save data as netcdf (not CF-conform!)
    print('\nSaving subset of gridded observations as netCDF...')
    grid.to_netcdf(dir_tmp+param+'_original_gridded_obs.nc')
    print(interp_method, distance, variogram_model, neighbor, type(neighbor))
    print('\nMerging station data with gridded observations...')
    print('Depending on the domain size this may take a while...')
    if param == 'pr':
        # call the function imported at the top and give the functions all the arguments needed.
        merge_obs = merge_no_grad(grid[param].data, param, st_data, lat, lon, st_lats, st_lons, dir_save_run,
                                  plot_opt=True, interp = interp_method, distance = distance,
                                  variogram_model = variogram_model, neighbor=neighbor)
    elif (param == 'tasmax') | (param == 'tasmin') | (param == 'rsds') | (param == 'sfcWind') | (param == 'hurs'):
        merge_obs = merge_hgt_grad(grid, param, df_st, lat, lon, st_lats, st_lons, st_height, topo.height.data, dir_save_run,
                                   plot_opt=True, interp = interp_method, distance = distance,
                                   variogram_model = variogram_model, neighbor=neighbor)

    # define a filename under which the new observations should be stored
    filename_obs = param+'_merged_observations'
    fn_obs_long = dir_save_run+filename_obs+'_1981-2010.nc'
    # save new observational dataset as netCDF (not CF-conform!)
    print('\nSaving new observational dataset as netCDF...')
    print(fn_obs_long)

    write_netcdf(merge_obs, param, lat, lon, 1981, 2010, dir_save_run, filename_obs)

# Get Climate Model Data, cut out the domain defined at the beginning and save data as netcdf
    print('\nReading historical model data...')
    if os.path.isfile(model_hist):
        mod_hist = xr.open_dataset(model_hist)
        # get modelname and calendar from metadata
        model_cal, model_name = get_ncattrs(model_hist)
    else:
        mod_hist = xr.open_mfdataset(model_hist+param+'*_historical_*.nc')
        model_cal, model_name = get_ncattrs(glob.glob(model_hist+param+'*_historical_*.nc')[0])
    mod_hist = cut_domain(mod_hist, lat_min, lat_max, lon_min, lon_max)

    ## model scenario data (rcp)
    print('\nReading model scenario data...')

    if start_y <= 2010:
        if end_y <= 2010:
            mod_rcp = mod_hist.copy(deep=True)
            if model_cal == '360_day':
                mod_rcp = mod_rcp.sel(time=slice(str(start_y)+'-01-01', str(end_y)+'-12-30'))
            else:
                mod_rcp = mod_rcp.sel(time=slice(str(start_y)+'-01-01', str(end_y)+'-12-31'))
        else:
            mod_rcp1 = mod_hist.copy(deep=True)
            if model_cal == '360_day':
                mod_rcp1 = mod_rcp1.sel(time=slice(str(start_y)+'-01-01', '2010-12-30'))
            else:
                mod_rcp1 = mod_rcp1.sel(time=slice(str(start_y)+'-01-01', '2010-12-31'))

            if os.path.isfile(model_rcp):
                mod_rcp2 = xr.open_dataset(model_rcp)
                model_cal, model_name = get_ncattrs(model_rcp)
            else:
                mod_rcp2 = xr.open_mfdataset(model_rcp+param+'*_rcp*_*.nc')
                model_cal, model_name = get_ncattrs(glob.glob(model_rcp+param+'*_rcp_*.nc')[0])

            if model_cal == '360_day':
                mod_rcp2 = mod_rcp2.sel(time=slice(str(start_y)+'-01-01', str(end_y)+'-12-30'))
            else:
                mod_rcp2 = mod_rcp2.sel(time=slice(str(start_y)+'-01-01', str(end_y)+'-12-31'))
            mod_rcp2 = cut_domain(mod_rcp2, lat_min, lat_max, lon_min, lon_max)

            mod_rcp = xr.concat([mod_rcp1, mod_rcp2], 'time')

    else:
        if os.path.isfile(model_rcp):
            mod_rcp = xr.open_dataset(model_rcp)
            model_cal, model_name = get_ncattrs(model_rcp)
        else:
            mod_rcp = xr.open_mfdataset(model_rcp+param+'*_rcp*_*.nc')
            model_cal, model_name = get_ncattrs(glob.glob(model_rcp+param+'*_rcp_*.nc')[0])

        if model_cal == '360_day':
            mod_rcp = mod_rcp.sel(time=slice(str(start_y)+'-01-01', str(end_y)+'-12-30'))
        else:
            mod_rcp = mod_rcp.sel(time=slice(str(start_y)+'-01-01', str(end_y)+'-12-31'))
        mod_rcp = cut_domain(mod_rcp, lat_min, lat_max, lon_min, lon_max)

    # check units in climate model
    if param == 'pr':
        if np.nanmean(mod_hist[param]) < 0.01:
            mod_hist[param] = mod_hist[param]*86400
        if np.nanmean(mod_rcp[param]) < 0.01:
            mod_rcp[param] = mod_rcp[param]*86400
    elif (param == 'tasmax') | (param == 'tasmin'):
        if np.nanmean(mod_hist[param]) > 250:
            mod_hist[param] = mod_hist[param] - 273.15
        if np.nanmean(mod_rcp[param]) > 250:
            mod_rcp[param] = mod_rcp[param] - 273.15

    # save cropped model data
    print('\nSaving subset of original historical model data as netCDF...')
    mod_hist_savename = dir_tmp+'/'+param+'_model_historical_subset_1981-2010.nc'
    print(mod_hist_savename)
    write_netcdf_model(mod_hist[param].data, param, lat, lon, 1981, 2010, mod_hist_savename, model_cal, model_name)

    print('\nSaving subset of original model data for period of bias correction as netCDF...')
    mod_rcp_savename = dir_tmp+'/'+param+'_model_bc_period_subset_'+str(start_y)+'-'+str(end_y)+'.nc'
    print(mod_rcp_savename)
    write_netcdf_model(mod_rcp[param].data, param, lat, lon, start_y, end_y, mod_rcp_savename, model_cal, model_name)


# Run Bias Correction
    print('\nRunning bias correction for every grid point in domain...')
    print('Depending on your domain size and grid resolution this will take a while...')
    print('The corrected grid point time series will be saved at '+dir_save_run+'BIASCORR_TS/')

    sdm_wrap(param, fn_obs_long, mod_hist_savename, mod_rcp_savename, dir_save_run, start_y, end_y)

# Reorder bias corrected data (grid point time series) to lat/lon field
# By default the final data will be stored in 10-year blocks.
    print('\nReordering and saving new bias corrected data to 3-dimensional array (time,lat,lon)...')

    reorder_data(param, dir_save_run+'BIASCORR_TS/', dir_save_run+'BIASCORR/', model_name, lat, lon, model_cal, start_y, end_y)

    ds = xr.open_mfdataset(dir_save_run+'BIASCORR/'+param+'*.nc')
    data_bc = ds[param]
    ds.close()

    ds = xr.open_dataset(mod_rcp_savename)
    data_bcp = ds[param]
    ds.close()




# plotting the results from biascorrection
    plot_result(param, data_bc.values, data_bcp.values, dir_save_run+'PLOTS/', lat, lon, st_lats, st_lons, start_y, end_y, 'biascorr','orig. model', 'improved biascorr.', False)
    print('\nYOU ARE DONE!')
    return
