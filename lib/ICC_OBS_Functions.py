# -*- coding: utf-8 -*-
"""
ICC-OBS Functions

Created on Thu Jun 14 15:44:48 2018
@author: mariaw
"""

import glob
import xarray as xr
import numpy as np
#import scipy as sp
#import scipy.interpolate
#import matplotlib.pyplot as plt
#import pykrige

# packages for writing netCDF
from netCDF4 import Dataset
from netCDF4 import date2num, num2date
from datetime import datetime, timedelta
import os.path


def get_ncattrs(fn_nc):
    nc_fid = Dataset(fn_nc, 'r')
    model_cal = str(nc_fid.variables['time'].calendar)
    model_name = str(nc_fid.getncattr('modelname'))
    return model_cal, model_name

def check_domain(topo, lat_min, lat_max, lon_min, lon_max):
    '''
    Check if specified domain is inside model domain.
    '''
    err = 0
    if (lat_min < topo.lat.min()) | (lat_max > topo.lat.max()) | (lon_min < topo.lon.min()) | (lon_max > topo.lon.max()):
        err = 1
    return err

def plausibility_check(param, st_data, st_lons, st_lats, lat_min, lat_max, lon_min, lon_max):
    '''
    check station values for plausibility (unrealistic values or undetected missing values in data)

    returns err = 0 if no wrong value is detected
    returns err != 0 and a warning message if a wrong value is detected
    '''
    err = 0

    if param == 'pr':
        if np.nanmin(st_data) < 0:
            err = 1
        if np.nanmax(st_data) > 500:
            err = 2

    elif (param == 'tasmax') | (param == 'tasmin'):
        if np.nanmin(st_data) < -50:
            err = 1
        if np.nanmax(st_data) > 50:
            err = 2

    elif param == 'rsds':
        if np.nanmin(st_data) < 0:
            err = 1
        if np.nanmax(st_data) > 450:
            err = 2

    elif (param == 'rh') | (param == 'hurs'):
        if np.nanmin(st_data) < 0:
            err = 1
        if np.nanmax(st_data) > 100:
            err = 2

    elif (param == 'ws10') | (param == 'sfcWind'):
        if np.nanmin(st_data) < 0:
            err = 1
        if np.nanmax(st_data) > 100:
            err = 2
    elif (lat_min < st_lats.min()) | (lat_max > st_lats.max()) | (lon_min < st_lons.min()) | (lon_max > st_lons.max()):
        err = 3
        print('ERROR: One or more stations are located outside the specified domain!')

    if err == 0:
        print ('Station data check ok!')
    elif err == 1:
        print ('ERROR: Unrealistic minimum value encounterd in station data!')
    elif err == 2:
        print ('ERROR: Unrealistic maximum value encounterd in station data!')

    return err

def cut_domain(ds, lat_min, lat_max, lon_min, lon_max):
    ds_cut = ds.where(((ds['lon'] >= lon_min) & (ds['lon'] <= lon_max) & (ds['lat'] >= lat_min) & (ds['lat'] <= lat_max)), drop=True)
    return ds_cut


def read_yearly_files(data_dir, fname, start_year, end_year, lat_min, lat_max, lon_min, lon_max, dy=1):
    # function reads multiple netCDF files into one dataset
    year = start_year
    while year < end_year:
    #for year in range(start_year, end_year+1):
        print(str(year))
        fil = glob.glob(data_dir+fname+'*'+str(year)+'*.nc')[0]
        ds = xr.open_dataset(fil)
        ds = cut_domain(ds, lat_min, lat_max, lon_min, lon_max)
        if year == start_year:
            param_all = ds
        else:
            param_all = xr.concat([param_all, ds], 'time')
        year = year+dy

    return param_all


def read_multiple_files(data_dir, fname, start_year, end_year, dyear, lat_min, lat_max, lon_min, lon_max):
    # function reads yearly netCDF files into one dataset

    for year in np.arange(start_year, end_year, dyear):
        print(str(year))
        fil = glob.glob(data_dir+fname+'*'+str(year)+'*.nc')[0]
        ds = xr.open_dataset(fil)
        ds = cut_domain(ds, lat_min, lat_max, lon_min, lon_max)
        if year == start_year:
            param_all = ds
        else:
            param_all = xr.concat([param_all, ds], 'time')

    return param_all


def write_netcdf(array3d, param_name, lat1d, lon1d, start_year, end_year, savedir, filename, cal='gregorian'):#, freq ='daily'):

#    fillval = -9999
#    array3d[np.isnan(array3d)] = fillval

    # create netCDF file
    savename = savedir+filename+'_'+str(start_year)+'-'+str(end_year)+'.nc'
    # check if file already exists - if yes, delete it before writing new file
    if os.path.exists(savename):
        os.remove(savename)

    # open new netCDF file in write mode:
    try:
        dataset = Dataset(savename,'w',format='NETCDF4_CLASSIC')
    except(IOError):
        os.remove(savename)
        dataset = Dataset(savename,'w',format='NETCDF4_CLASSIC')
    # create dimensions

    dataset.createDimension("time",None)
    dataset.createDimension("lat", len(lat1d))
    dataset.createDimension("lon",len(lon1d))

    # create variables
    times = dataset.createVariable("time","f8",("time",))
    lats = dataset.createVariable("lat","f4",("lat",))
    lons = dataset.createVariable("lon","f4",("lon",))
    var = dataset.createVariable(param_name, "f4", ("time","lat","lon",), fill_value = -9999)
    crs = dataset.createVariable('crs', 'i', ())
    #prec = dataset.createVariable('pr', "f4", ("time","y","x",))

    # add attributes
    times.units = 'days since 1950-01-01T00:00:00Z'
    times.calendar = cal
    times.long_name = 'time'
    times.axis = 'T'
    times.standard_name = 'time'

    lats.units = 'degrees_north'
    lats.long_name = 'latitude'
    lats.standard_name = 'latitude'

    lons.units = 'degrees_east'
    lons.long_name = 'longitude'
    lons.standard_name = 'longitude'

    if (param_name == 'pr') | (param_name == 'rr'):
        var.units = 'mm'
        var.long_name = 'total daily precipitation'
        var.standard_name = 'precipitation_amount'
        dataset.title = 'daily precipitation amount'
        print('...writing pr')

    elif (param_name == 'tasmax') | (param_name == 'tmax'):
        var.units = 'degree_Celsius'
        var.long_name = 'daily maximum near-surface air temperature'
        var.standard_name = 'air_temperature'
        dataset.title = 'daily maximum temperature'
        print('...writing tmax')

    elif (param_name == 'tasmin') | (param_name == 'tmin'):
        var.units = 'degres_Celsius'
        var.long_name = 'daily minimum near-surface air temperature'
        var.standard_name = 'air_temperature'
        dataset.title = 'daily minimum temperature'
        #dataset.comment = 'Merged gridded observations of daily precipitation amount using DANUBECLIM as primary source and E-OBS (Version 16.0) data (regridded with ESMF_RegridWeightGen) as secondary source.'
        print('...writing tmin')

    elif (param_name == 'rsds'):
        var.units = 'W m-2'
        var.long_name = 'surface downwelling shortwave flux'
        var.standard_name = 'surface_downwelling_shortwave_flux_in_air'
        dataset.title = 'daily mean global radiation'
        print('...writing rsds')

    elif (param_name == 'sfcWind') | (param_name == 'ws'):
        var.units = 'm s-1'
        var.long_name = 'daily mean 10-m wind speed'
        var.standard_name = 'wind_speed'
        dataset.title='daily mean 10-m wind speed'

    elif (param_name == 'hurs') | (param_name == 'rh'):
        var.units = 'percent'
        var.long_name = 'daily mean relative humidity'
        var.standard_name = 'relative_humidity'
        dataset.title='daily mean relative humidity'

    else:
        print ('unknown parameter - edit information in function!')
        print ('aborting... no netCDF file saved!')
        dataset.close()
        os.remove(savename)
        return

    var.grid_mapping = 'latitude_longitude'

    crs.grid_mapping_name = "latitude_longitude"
    crs.longitude_of_prime_meridian = 0.0
    crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257223563
    crs.comment = 'Latitude and longitude on the WGS 1984 datum'


    # write data to netCDF variable
    array3d[np.isnan(array3d)] = -9999
    var[:] = array3d
    lats[:] = lat1d
    lons[:] = lon1d

    if str(cal) == '365_day':
        d = np.arange((start_year-1950)*365, (end_year-1950+1)*365)
        dates = num2date(d, 'days since 1950-01-01T00:00:00Z', calendar=cal)
    elif str(cal) == '360_day':
        d = np.arange((start_year-1950)*360, (end_year-1950+1)*360)
        dates = num2date(d, 'days since 1950-01-01T00:00:00Z', calendar=cal)
    else:
        dates = [datetime(start_year,1,1)+k*timedelta(days=1) for k in range(array3d.shape[0])]
    times[:] = date2num(dates, units=times.units, calendar=times.calendar)

    # global attributes
    dataset.project="Climaproof, funded by the Austrian Development Agency (ADA) and co-funded by the United Nations Environmental Programme (UNEP)"
    dataset.source = "Observational dataset created with the ICC-OBS Tool (Institute of Meteorology, University of Natural Resources and Life Sciences, Vienna, Austria)"
    dataset.conventions = "CF-1.6"

    dataset.close()

    print('Writing your data to netCDF was successful!')

    return

def write_netcdf_model(array3d, param_name, lat1d, lon1d, start_year, end_year, savename, cal, model_name):#, freq ='daily'):

#    fillval = -9999
#    array3d[np.isnan(array3d)] = fillval

    # create netCDF file
    # check if file already exists - if yes, delete it before writing new file
    if os.path.exists(savename):
        os.remove(savename)

    # open new netCDF file in write mode:
    try:
        dataset = Dataset(savename,'w',format='NETCDF4_CLASSIC')
    except(IOError):
        os.remove(savename)
        dataset = Dataset(savename,'w',format='NETCDF4_CLASSIC')
    # create dimensions

    dataset.createDimension("time",None)
    dataset.createDimension("lat", len(lat1d))
    dataset.createDimension("lon",len(lon1d))

    # create variables
    times = dataset.createVariable("time","f8",("time",))
    lats = dataset.createVariable("lat","f4",("lat",))
    lons = dataset.createVariable("lon","f4",("lon",))
    var = dataset.createVariable(param_name, "f4", ("time","lat","lon",), fill_value = -9999)
    crs = dataset.createVariable('crs', 'i', ())
    #prec = dataset.createVariable('pr', "f4", ("time","y","x",))

    # add attributes
    times.units = 'days since 1950-01-01T00:00:00Z'
    times.calendar = cal
    times.long_name = 'time'
    times.axis = 'T'
    times.standard_name = 'time'

    lats.units = 'degrees_north'
    lats.long_name = 'latitude'
    lats.standard_name = 'latitude'

    lons.units = 'degrees_east'
    lons.long_name = 'longitude'
    lons.standard_name = 'longitude'

    if (param_name == 'pr'):
        var.units = 'mm'
        var.long_name = 'total daily precipitation'
        var.standard_name = 'precipitation_amount'
        dataset.title = 'daily precipitation amount'
        print('...writing pr')

    elif (param_name == 'tasmax'):
        var.units = 'degree_Celsius'
        var.long_name = 'daily maximum near-surface air temperature'
        var.standard_name = 'air_temperature'
        dataset.title = 'daily maximum temperature'
        print('...writing tmax')

    elif (param_name == 'tasmin'):
        var.units = 'degres_Celsius'
        var.long_name = 'daily minimum near-surface air temperature'
        var.standard_name = 'air_temperature'
        dataset.title = 'daily minimum temperature'
        print('...writing tmin')

    elif (param_name == 'rsds'):
        var.units = 'W m-2'
        var.long_name = 'surface downwelling shortwave flux'
        var.standard_name = 'surface_downwelling_shortwave_flux_in_air'
        dataset.title = 'daily mean global radiation'
        print('...writing rsds')

    elif (param_name == 'sfcWind'):
        var.units = 'm s-1'
        var.long_name = 'daily mean 10-m wind speed'
        var.standard_name = 'wind_speed'
        dataset.title='daily mean 10-m wind speed'

    elif (param_name == 'hurs'):
        var.units = 'percent'
        var.long_name = 'daily mean relative humidity'
        var.standard_name = 'relative_humidity'
        dataset.title='daily mean relative humidity'

    else:
        print ('unknown parameter - edit information in function!')
        print ('aborting... no netCDF file saved!')
        dataset.close()
        os.remove(savename)
        return

    var.grid_mapping = 'latitude_longitude'

    crs.grid_mapping_name = "latitude_longitude"
    crs.longitude_of_prime_meridian = 0.0
    crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257223563
    crs.comment = 'Latitude and longitude on the WGS 1984 datum'


    # write data to netCDF variable
    array3d[np.isnan(array3d)] = -9999
    var[:] = array3d
    lats[:] = lat1d
    lons[:] = lon1d

    if str(cal) == '365_day':
        d = np.arange((start_year-1950)*365, (end_year-1950+1)*365)
        dates = num2date(d, 'days since 1950-01-01T00:00:00Z', calendar=cal)
    elif str(cal) == '360_day':
        d = np.arange((start_year-1950)*360, (end_year-1950+1)*360)
        dates = num2date(d, 'days since 1950-01-01T00:00:00Z', calendar=cal)
    else:
        dates = [datetime(start_year,1,1)+k*timedelta(days=1) for k in range(array3d.shape[0])]
    times[:] = date2num(dates, units=times.units, calendar=times.calendar)

    # global attributes
    dataset.project="Climaproof, funded by the Austrian Development Agency (ADA) and co-funded by the United Nations Environmental Programme (UNEP)"
    dataset.source = "Climate model subset created with the ICC-OBS Tool (Institute of Meteorology, University of Natural Resources and Life Sciences, Vienna, Austria)"
    dataset.model = model_name
    dataset.conventions = "CF-1.6"

    dataset.close()

    print ('Writing your data to netCDF was successful!')

    return

def plot_result(param, data_bc, data_hist, plot_dir, lat, lon, st_lats, st_lons, start_y, end_y, plotname, title1, title2, scatter):

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    from mpl_toolkits.axes_grid1 import AxesGrid
    import cartopy as cpy
    import cartopy.crs as ccrs
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    years = int(end_y) - int(start_y) + 1

    borders = cpy.io.shapereader.Reader("lib/ne_50m_admin_0_countries.shp")

    if not os.path.exists(plot_dir): # creates directory if not existing
        os.makedirs(plot_dir)

    if (param == 'tasmax') | (param == 'tasmin'):
        data_hist = data_hist.mean(axis=0)
        data_bc = data_bc.mean(axis=0)
        cmap = 'YlOrRd'

    elif (param == 'rsds') | (param == 'hurs'):
        data_hist = data_hist.mean(axis=0)
        data_bc = data_bc.mean(axis=0)
        cmap = 'viridis'

    elif param == 'sfcWind':
        data_hist = data_hist.mean(axis=0)
        data_bc = data_bc.mean(axis=0)
        cmap = 'RdYlGn_r'

    if param != 'pr':
        vmin = np.floor(np.nanmin([np.nanmin(data_bc), np.nanmin(data_hist)]))
        vmax = np.ceil(np.nanmax([np.nanmax(data_bc), np.nanmax(data_hist)]))
        cmap_diff = 'bwr'

        if (vmax - vmin) <= 5:
            step = 0.5
        elif (vmax - vmin) <= 10:
            step = 1.
        elif (vmax - vmin) <= 20:
            step = 2.
        elif (vmax - vmin) <= 50:
            step = 5.
        elif (vmax - vmin) > 50:
            step = 10.

        levels = np.arange(vmin, vmax+step, step)

    if param == 'pr':
        cmap_diff = 'bwr_r'
        cmap = 'YlGnBu'
        data_hist = data_hist.sum(axis=0)/years
        data_bc = data_bc.sum(axis=0)/years
        vmin = np.floor(np.nanmin([np.nanmin(data_bc), np.nanmin(data_hist)]))
        vmax = np.ceil(np.nanmax([np.nanmax(data_bc), np.nanmax(data_hist)]))
        vmin = np.floor(vmin / 100) * 100
        vmax = np.ceil(vmax / 100) * 100

        if (vmax - vmin) <= 200:
            step = 30
        elif (vmax - vmin) <= 400:
            step = 60
        elif (vmax - vmin) <= 800:
            step = 100
        elif (vmax - vmin) <= 1600:
            step = 200
        elif (vmax - vmin) <= 2000:
            step = 250
        elif (vmax - vmin) <= 3000:
            step = 350
        else:
            step = 500

        levels = np.arange(vmin, vmax+step, step)

    projection = ccrs.PlateCarree()

    axes_class = (GeoAxes,
                  dict(map_projection=projection))

    delta_lat = np.ceil((lat.max() - lat.min()) * 1.2)
    delta_lon = np.ceil((lon.max() - lon.min()) * 3.5)
    print(delta_lat, delta_lon)

    fig = plt.figure(1, figsize=(delta_lon,delta_lat))
    fig.suptitle(param+' ('+str(start_y)+' - '+str(end_y)+')', fontsize=14)
    ax = AxesGrid(fig, 111, axes_class=axes_class,
                    nrows_ncols=(1, 3),
                    axes_pad=0.6,
                    share_all=True,
                    cbar_location='right',
                    cbar_mode='each',
                    cbar_pad=0.1,
                    cbar_size='3%',
                    label_mode='')  # note the empty label_mode

    for i in ax:
        for land in borders.geometries():
            i.add_geometries([land], crs=projection, facecolor='none', linewidth=0.7, edgecolor="black")

    if lat.all() > 0:
        lat_hem = 'N'
    elif lat < 0:
        lat_hem = 'S'
    else:
        lat_hem = ''

    if lon.all() > 0:
        lon_hem = 'E'
    elif lon < 0:
        lon_hem = 'W'
    else:
        lon_hem = ''

    degree = '$^\circ$'

    if lat.max() - lat.min() <= 5:
        lat_num = 5
    elif lat.max() - lat.min() <= 7:
        lat_num = 7
    elif lat.max() - lat.min() <= 9:
        lat_num = 9
    else: lat_num = 11

    if lon.max() - lon.min() <= 5:
        lon_num = 5
    elif lon.max() - lon.min() <= 7:
        lon_num = 7
    elif lon.max() - lon.min() <= 9:
        lon_num = 9
    else: lon_num = 11


    gl = ax[0].gridlines(crs=projection, draw_labels=True,
                  linewidth=0.5, color='k', alpha=0.9, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator(np.linspace(lon.min(), lon.max(), num=lon_num, endpoint=True))
    gl.ylocator = mticker.FixedLocator(np.linspace(lat.min(), lat.max(), num=lat_num, endpoint=True))
    gl.xformatter = mticker.FormatStrFormatter("%.1f"+degree+lon_hem)
    gl.yformatter = mticker.FormatStrFormatter("%.1f"+degree+lat_hem)
    gl.ylabel_style = {'size': 7, 'color': 'black'}
    gl.xlabel_style = {'size': 7, 'color': 'black', 'rotation': 30}

    gl = ax[1].gridlines(crs=projection, draw_labels=True,
                  linewidth=0.5, color='k', alpha=0.9, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator(np.linspace(lon.min(), lon.max(), num=lon_num, endpoint=True))
    gl.ylocator = mticker.FixedLocator(np.linspace(lat.min(), lat.max(), num=lat_num, endpoint=True))
    gl.xformatter = mticker.FormatStrFormatter("%.1f"+degree+lon_hem)
    gl.xlabel_style = {'size': 7, 'color': 'black', 'rotation': 30}

    gl = ax[2].gridlines(crs=projection, draw_labels=True,
                  linewidth=0.5, color='grey', alpha=0.9, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator(np.linspace(lon.min(), lon.max(), num=lon_num, endpoint=True))
    gl.ylocator = mticker.FixedLocator(np.linspace(lat.min(), lat.max(), num=lat_num, endpoint=True))
    gl.xformatter = mticker.FormatStrFormatter("%.1f"+degree+lon_hem)
    gl.xlabel_style = {'size': 7, 'color': 'black', 'rotation': 30}

    ax[0].set_extent((lon.min(), lon.max(), lat.min(), lat.max()), crs=projection)
    f1 = ax[0].contourf(lon,lat,data_hist, transform=projection, levels=levels, cmap=cmap)
    if scatter == True:
        ax[0].scatter(st_lons,st_lats, s=20, transform=projection, c='k', marker='^')
    ax[0].set_title(title1)
    ax.cbar_axes[0].colorbar(f1, ticks=levels)

    ax[1].set_extent((lon.min(), lon.max(), lat.min(), lat.max()), crs=projection)
    f2 = ax[1].contourf(lon,lat,data_bc, transform=projection, levels=levels, cmap=cmap)
    if scatter == True:
        ax[1].scatter(st_lons,st_lats, s=20, transform=projection, c='k', marker='^')
    ax[1].set_title(title2)
    ax.cbar_axes[1].colorbar(f2, ticks=levels)

    vmin = np.floor(np.nanmin(data_bc - data_hist))
    vmax = np.ceil(np.nanmax(data_bc - data_hist))

    if np.abs(vmin) > np.abs(vmax):
        vmax = np.abs(vmin)
    else:
        vmin = - vmax

    if (param == 'pr'):
        if (vmax <=100):
            step = 50
        elif vmax <= 200:
            step = 100
        elif vmax <= 1000:
            step = 250
        elif vmax <= 2000:
            step = 500
        elif vmax <= 3000:
            step = 1000

    levels = np.linspace(vmin, vmax, num=11, endpoint=True)
    f3 = ax[2].contourf(lon, lat,(data_bc - data_hist), transform=projection, levels=levels, cmap=cmap_diff)
    ax[2].set_extent((lon.min(), lon.max(), lat.min(), lat.max()), crs=projection)
    if scatter == True:
        ax[2].scatter(st_lons,st_lats, s=20, transform=projection, c='k', marker='^')
    ax[2].set_title('difference\n'+title2 +' - '+ title1)
    ax.cbar_axes[2].colorbar(f3, ticks=levels)

    plt.subplots_adjust(top=0.95,
                        bottom=0.05,
                        left=0.05,
                        right=0.92,
                        hspace=0.2,
                        wspace=0.15)
    plt.savefig(plot_dir+param+'_'+plotname+'.png', dpi=400)
    plt.close()
    return
