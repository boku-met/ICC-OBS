# -*- coding: utf-8 -*-

#import xarray as xr
import numpy as np
import scipy as sp
#import scipy.interpolate
import matplotlib.pyplot as plt
import pykrige
import os
from ICC_OBS_Functions import plot_result
from metpy.interpolate import inverse_distance_to_grid
from metpy.interpolate import remove_nan_observations




#%% KRIGING INTERPOLATION
def kriging(st_lons, st_lats, res_data, lon, lat, var_model=None, nlags=10, var_dict = None):
    '''
    This function calls the kriging routine to interpolate the residuals
    'res_data', that are defined at the points 'st_lons' and 'st_lats' to
    a grid defined by 'lat' and 'lon'.

    https://pykrige.readthedocs.io/en/latest/generated/pykrige.ok.OrdinaryKriging.htm

    var_model:  variogramm model (gaussian, sperical, exponential or linear)

    nlags:      number of bins, the semivariogramm is divided into (higher = more detailed but slower)

    var_dict:    dictionary defining the parameters of the semivariogramm
                vardict = {'sill': s, 'range': r, 'nugget': n} or None (automatic calculation of parameters)
    '''
    # define variogram_model for kriging
    #var_model = 'spherical'

    OK_res = pykrige.ok.OrdinaryKriging(st_lons, st_lats, res_data, variogram_model=var_model, 
                                        variogram_parameters=var_dict, weight=True, verbose=False, 
                                        enable_plotting=False, nlags = nlags, coordinates_type='geographic')
    res_int, _ = OK_res.execute('grid', lon, lat)
#
#    UK_res = pykrige.uk.UniversalKriging(st_lons, st_lats, res_data, variogram_model=var_model,
#                      drift_terms=['regional_linear'])
#    res_int_uk, _ = UK_res.execute('grid', lon, lat)
    return res_int

#%% INVERSE DISTANCE INTERPOLATION
def haversine_dist(lat1,lat2,lon1,lon2):
    '''
    This function calculates the distance between lat/lon points
    '''
    # Distance between Grid of Lat and Lons and Location
    # lat_grid, lat_test, lon_grid, lon_test
    lat1 = sp.array(lat1)
    lon1 = sp.array(lon1)

    r=6378.000
    rad=lambda grad: grad*sp.pi/180.
    d=2*r*sp.arcsin(sp.sqrt(sp.sin((rad(lat2)-rad(lat1))/2.)**2+sp.cos(rad(lat1))*sp.cos(rad(lat2))*sp.sin((rad(lon2)-rad(lon1))/2.)**2))
    return d

#def calc_weights(st_lats, st_lons, grid_lats, grid_lons, weights, epsilon, low_lim):
#    '''
#    Calculate weight matrix for inverse distance interpolation
#
#    Returns:
#    weight matrix (w) and distance matrix (dist)
#
#    *************************************************
#    INPUT:
#
#    st_lats, st_lons        1d arrays with latitude and longitude values of stations
#
#    grid_lats, grid_lons    1d arrays with latitude and longitude values of grid
#
#    weights                 function defining the decay of the weights with distance
#                            available functions: gaussian, idw, linear or quadratic
#
#    epsilon                 parameter defining the slope of the weight function
#
#    low_lim                 threshold below which the weights are set to 0   
#    '''   
#    epsilon = float(epsilon)
#    nst = len(st_lats) # number of stations
#
#    w = np.zeros((nst, len(grid_lats)*len(grid_lons)))
#    dist = np.zeros((nst, len(grid_lats)*len(grid_lons)))
#    for i in range(nst):
#        ij = 0
#        for j in range(len(grid_lats)):
#            for k in range(len(grid_lons)):
#                #ij = (j+1)*(k+1)-1
#                dist[i, ij] = haversine_dist(st_lats[i], grid_lats[j], st_lons[i], grid_lons[k])
#                if weights == 'idw':
#                    w[i, ij] = (np.nanmax([0., 100. - dist[i, ij]]) / (100.*dist[i, ij]))**2
##                    w[i, ij] = 1.0/np.sqrt((dist[i, ij]/epsilon)**2 + 1)
##                    if w[i,ij] <= low_lim:
##                        w[i, ij] = 0
#                elif weights == 'gaussian':
#                    w[i, ij] = np.exp(-(dist[i, ij]/epsilon)**2)
#                    if w[i,ij] <= low_lim:
#                        w[i, ij] = 0
#                elif weights == 'linear':
#                    w[i, ij] = 1.0-(epsilon*dist[i, ij])
#                    if w[i,ij] <= 0:
#                        w[i, ij] = 0
#                elif weights == 'quadratic':
#                    w[i, ij] = 1.0-(epsilon*dist[i, ij]**2)
#
#                ij = ij+1
#
#    return w, dist

#def idw_interpolation(data, st_lats, st_lons, grid_lats, grid_lons, weights=None, distance = None, low_lim = None):
#    '''
#    Inverse Distance interpolation of station data to grid
#
#    Returns:
#    3d array with interpolated station data on grid (dimensions: time, lat, lon)
#
#    ******************************************************
#    INPUT:
#
#    data                    time series of data at stations (residuals) - dimensions: time, # of stations
#
#    st_lats, st_lons        1d arrays with latitude and longitude values of stations
#
#    grid_lats, grid_lons    1d arrays with latitude and longitude values of grid
#
#    weights                 function defining the decay of the weights with distance
#                            gaussian (default), idw, linear or quadratic
#
#    distance                distance (in km) at which the weights should be 0 (default: 50km)
#
#    low_lim                 threshold below which the weights are set to 0 (default: 0.05)   
#        
#    '''
#    # calculate epsilon (for weight function) according to distance
#    epsilon = get_epsilon(distance, weights, low_lim)
#
#    # calculate weight matrix
#    w, dist = calc_weights(st_lats, st_lons, grid_lats, grid_lons, weights = weights, epsilon = epsilon, low_lim = low_lim)
#
#    # interpolation (loop over every day)
#    grid_interp = np.zeros((np.shape(data)[0], len(grid_lats), len(grid_lons)))
#    for t in range(np.shape(data)[0]):
#        w_tmp = np.copy(w)
#        arr = data[t,:]
#        # check for nans in station values
#        nan_idx = np.argwhere(np.isnan(arr))
#        arr[np.isnan(arr)] = 0
#        if len(nan_idx) > 0:
#            arr[nan_idx] = 0
#            w_tmp[nan_idx,:] = 0
#
#        grid_vals = np.dot(arr, w_tmp)
#        
#        # rearange vector to 2d array
#        grid_interp[t,:,:] = np.reshape(grid_vals, (len(grid_lats), len(grid_lons)))
#
#    return grid_interp

#def get_epsilon(dist, func, low_lim):
#    if func == 'idw':
#        epsilon = dist/np.sqrt(1./low_lim**2 -1.)
#    elif func == 'gaussian':
#        epsilon = dist/np.sqrt(np.abs(np.log(low_lim)))
#    elif func == 'quadratic':
#        epsilon = 1./(dist**2)
#    elif func == 'linear':
#        epsilon = 1./dist
#    return epsilon

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


#%% LINEAR REGRESSION

def linreg(data, topo, plot_opt=False):
    # remove nans
    nans = np.isnan(data)
    data_reg = data[~nans]
    topo_reg = topo[~nans]
    # linear regression
    A = np.vstack([topo_reg, np.ones(len(data_reg))]).T
    linreg = np.linalg.lstsq(A, data_reg, rcond=None)
    gradient ,constant = linreg[0]

    if plot_opt == True:
        #testplot of regression
        plt.figure()
        plt.scatter(topo, data)
        plt.plot(gradient*np.arange(0,2200,1)+constant, color='red')

    return gradient , constant


#%%  
def merge_no_grad(grid_data, param, st_data, lat, lon, st_lats, st_lons, savedir, plot_opt = False, interp = None, distance = None, variogram_model = None, neighbor=None):

    '''
    This function merges gridded observations with additional stationdata
    with the following steps:

    1.  Interpolation of the gridded observations to the stations using nearest
        neighbour interpolation method

    2.  Calculate residuals (difference between station data and the gridded
        observations interpolated to the stations) at the stations

    3.  Interpolate the residuals to the grid defined by 'lat' and 'lon' (grid
        of the gridded observations)
        Interpolation algorithm: Ordinary kriging
        https://pykrige.readthedocs.io/en/latest/generated/pykrige.ok.OrdinaryKriging.html

    4.  Merge the interpolated residuals with the original gridded observations
        to get the new gridded observations

    5.  Optional: Plot mean original and new observations for comparison

    ****************************************************
    Input data:

    grid_data   3d array of gridded observations; dimensions: (time,lat,lon)

    param       string with name of the parameter

    st_data     2d array with the station timeseries; dimensions: (time, station)

    lat         1d array with the latitude values of grid_data

    lon         1d array with the latitude values of grid_data

    st_lats     1d array with the latitude values of the stations

    st_lons     1d array with the longitude values of the stations

    plot_opt    default: False; set True to create simple plots of mean original and
                new observations

    interp      interpolation algorithm: kriging (default) or inverse distance (default
                if number of stations <=10)

    Parameters for inverse distance interpolation:
    
    neighbor        Minimum number of neighbors needed to perform barnes or 
                    cressman interpolation for a point. Default is 3.
                    See also: https://unidata.github.io/MetPy/latest/api/generated/metpy.interpolate.inverse_distance_to_grid.html

    distance        distance (in km) at which the weights should be 0 (default: 50km)
    

    Parameters for kriging:
        
    variogram_model         variogram model for the kriging interpolation:
                            gaussian(default), spherical, exponential or linear
     
    Output data:

    grid_merge  3d array with the merged observations; dimensions: (time,lat,lon)
    '''

    nst = np.shape(st_data)[1]
    # 1. Interpolate gridded observations to stations (nearest neighbour interpolation)
    grid_at_st = np.zeros((np.shape(grid_data)[0], nst))
    # loop over every day
    for i in range(np.shape(grid_data)[0]):
        grid_at_st[i, :] = sp.interpolate.interpn((lat, lon), grid_data[i, :, :], (st_lats, st_lons), method='linear')

    # 2. Calculate residuals at stations: station data - gridded observations interpolated to stations
    res_data = st_data - grid_at_st

    # 3. Interpolate residuals (ordinary kriging to 0.01° grid or inverse distance)

    res_int = np.zeros_like(grid_data)
    lons, lats = np.meshgrid(lon, lat)
        
    if (interp == 'idw'):
        #res_int = idw_interpolation(res_data, st_lats, st_lons, lat, lon, weights = weights, distance = distance, low_lim = low_lim)
        for x, val in enumerate(res_data):
            #print(val, x)
            lon_masked, lat_masked, val_masked = remove_nan_observations(st_lons, st_lats, val)
            res_int[x,:,:] = inverse_distance_to_grid(lon_masked, lat_masked, val_masked, lons, lats, r=distance, min_neighbors=neighbor )
    elif interp == 'kriging':
        # loop over evey day
        for i in range(np.shape(grid_data)[0]):
            # check for nans in observations
            mask = ~np.isnan(res_data[i, :])
            # check if all observations on this day are zero
            if np.all(res_data[i, mask] == 0.):
                res_int[i, :, :] = 0.
               # print 'all zeros at day ', str(i)
                continue
            else:
                res_int[i,:,:] = kriging(st_lons[mask], st_lats[mask], res_data[i, mask], lon, lat, var_model = variogram_model)
    
    res_int[np.isnan(res_int)] = 0    
    # final merged grid = original observation grid + interpolated residuals
    grid_merge = grid_data + res_int
    # check for values < 0 and set them to zero
    grid_merge[grid_merge < 0] = 0

    if plot_opt == True:
        plot_result(param, grid_merge, grid_data, savedir+'/PLOTS/', lat, lon, st_lats, st_lons, '1981', '2010', interp,'original obs.', 'improved obs.', scatter=True)

    return grid_merge

#%%
def merge_hgt_grad(grid, param, st_data, lat, lon, st_lats, st_lons, st_height, topo_height, savedir, plot_opt=False, interp = None, distance = None, variogram_model = None, neighbor=None):
    
    '''
    This function merges gridded temperature observations with additional
    stationdata with the following steps:
 
    1.  Calculate monthly height gradient of gridded data with
        linear regression

    2.  Remove height dependency from gridded and station data

    3.  Interpolation of the gridded observations to the stations using nearest
        neighbour interpolation method

    4.  Calculate residuals (difference between station data and the gridded
        observations interpolated to the stations) at the stations

    5.  Interpolate the residuals to the grid defined by 'lat' and 'lon' (grid
        of the gridded observations)

        Default interpolation algorithm:
            Inverse Distance Weighting Interpolation (idw)
        Other method:
            Ordinary Kriging
            https://pykrige.readthedocs.io/en/latest/generated/pykrige.ok.OrdinaryKriging.html

    6.  Merge interpolated residuals with grid (withoud height dependency)
    
    7.  Add height dependency back to the merged grid
        
    8.  Optional: Plot mean original and new observations for comparison

    ****************************************************
    Input data:

    grid        dataset of gridded observations; dimensions: (time,lat,lon)

    param       string with name of the parameter (tasmax or tasmin)

    st_data     dataframe with the station timeseries; dimensions: (time, station)

    lat         1d array with the latitude values of grid_data

    lon         1d array with the latitude values of grid_data

    st_lats     1d array with the latitude values of the stations

    st_lons     1d array with the longitude values of the stations

    st_height   1d array with the height values of the stations

    topo_height 2d array with gridded height data; dimensions: (lat, lon)

    plot_opt    default: False; set True to create simple plots of mean original and
                new observations

    interp      specifies the interpolation algorithm that should be used to
                interpolate the residuals
                default: 'idw'; other options: 'kriging', 'rbf'

    Parameters for inverse distance interpolation:

    neighbor        Minimum number of neighbors needed to perform barnes or 
                    cressman interpolation for a point. Default is 3.
                    See also: https://unidata.github.io/MetPy/latest/api/generated/metpy.interpolate.inverse_distance_to_grid.html

    distance        distance (in km) at which the weights should be 0 (default: 50km)

    Parameters for kriging:
        
    variogram_model         variogram model for the kriging interpolation:
                            gaussian(default), spherical, exponential or linear
    
    Output data:

    grid_merge  3d array with the merged observations; dimensions: (time,lat,lon)

    '''

#%%
    # 1. Remove height dependency for calculating residuals
    # 1.1. Get mean height gradient from gridded and station data

    # Calculate mean monthly height gradient and detrend data

    # add coordinate month to dataframe (for monthly selection of data)
    grid.coords['month'] = grid['time.month']

    # empty array that gets filled with monthly mean gradient for every gridpoint
    grad = np.zeros((12,1))*np.nan
    c = np.zeros((12,1))*np.nan

    grid_det = grid.copy(deep=True)
    topo_1d = topo_height.reshape(np.size(topo_height),1)

    st_det = st_data.copy(deep=True)

    for month in range(0,12):

        # get data for month
        month_data = grid[param][grid['month']==month+1,:,:]
        month_mean = month_data.mean(dim='time')
        grid_1d = month_mean.data.reshape(np.size(month_mean),1)

        grad[month], c[month] = linreg(grid_1d, topo_1d, plot_opt=False)

        # station data
        month_st = st_data[st_data.index.month==month+1]

        # 2. Remove height dependency from data
        grid_det[param][grid['month']==month+1,:,:] = month_data - (grad[month]* topo_height)
        st_det[st_data.index.month==month+1] = month_st - (grad[month]* st_height)
    
#%%    
    # 3. Interpolate gridded observations to stations (nearest neighbour interpolation)
    grid_at_st = np.zeros((grid[param].shape[0], st_data.shape[1]))
    # loop over every day
    for ii in range(grid[param].shape[0]):
        grid_at_st[ii, :] = sp.interpolate.interpn((lat, lon), grid_det[param][ii, :, :].data, (st_lats, st_lons), method='linear')

    # 4. Calculate residuals at stations: station data - gridded observations interpolated to stations
    res_data = st_det.values - grid_at_st

 
#%%    
    # 5. Interpolate residuals
    res_int = np.zeros(grid[param].shape) #np.zeros_like(grid_data)
    lons, lats = np.meshgrid(lon, lat)
    
    if (interp == 'idw'):
        # inverse distance interpolation (for small number of stations)
        for x, val in enumerate(res_data):
            # print(val, x)
            lon_masked, lat_masked, val_masked = remove_nan_observations(st_lons, st_lats, val)
            res_int[x,:,:] = inverse_distance_to_grid(lon_masked, lat_masked, val_masked, lons, lats, r=distance, min_neighbors=neighbor )

    else:
        for i in range(grid[param].shape[0]):
            # check for nans in data (missing values)
            mask = ~np.isnan(res_data[i, :])
            if interp == 'rbf':
                intfunc_rbf = sp.interpolate.Rbf(st_lons[mask], st_lats[mask], res_data[i,mask], function='gaussian')
                res_int[i,:,:] = intfunc_rbf(lons, lats)

            elif interp == 'kriging':
                # ordinary kriging
                res_int[i,:,:] = kriging(st_lons[mask], st_lats[mask], res_data[i, mask], lon, lat, var_model = variogram_model)           
    
            else:
                print('interpolation algorithm '+str(interp)+' unknown!')
                print('... aborting interpolation - no merged observations created!')
                return np.nan

    # 6. Calculate final field

    res_int[np.isnan(res_int)] = 0
    # 6.1 Merge residuals with gridded data
    grid_merge_det = grid_det[param].data + res_int
    
    # 6.2 Add height dependency
    grid_merge = np.ones_like(grid_merge_det)*np.nan
    
    for month in range(0,12):
        grid_merge[grid['month']==month+1,:,:] = grid_merge_det[grid['month']==month+1,:,:] + (grad[month]*topo_height)
        
    # check for unrealistic values
    if (param == 'rsds') | (param == 'sfcWind') | (param == 'hurs'):
        grid_merge[grid_merge < 0] = 0

    if param == 'hurs':
        grid_merge[grid_merge < 10] = 10
        grid_merge[grid_merge > 100] = 100
    print(grid_merge)
    if plot_opt == True:
        plot_result(param, grid_merge, grid[param].values, savedir+'/PLOTS/', lat, lon, st_lats, st_lons, '1981', '2010', interp,'original obs.', 'improved obs.', scatter=True)
    
    return grid_merge
