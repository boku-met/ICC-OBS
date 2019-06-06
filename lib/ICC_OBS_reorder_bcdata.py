# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
from netCDF4 import Dataset
from netCDF4 import date2num, num2date
from datetime import date, datetime, timedelta
import os.path
import glob

def reorder_data(param, dir_ts, dir_save, model_name, lat, lon, cal, startyear, maxyear):
    """    
    This function reorders the bias corrected time series of each grid point to
    2-dimensional lat-lon fields (same format as original model data).
    
    ****************
    INPUT:    
    
    param:      name of parameter (string)
    dir_ts:     directory with bias corrected grid point time series   
    dir_save:   directory where to save the data
    model_name: name of the climate model - for saving the data
    lat, lon:   1d arrays of latitude and longitude values of domain
    cal:        calendar type of model data as it is defined in the original 
                model data (e.g. 'proleptic_gregorian')
    startyear:  start year of grid point time series 
    maxyear:    end year of grid point time series
    """  
    
    # define block of years to save (10 = files with 10-years of data)
    dy = 10

    lat1d = lat
    lon1d = lon       
    year = startyear 

    if not os.path.exists(dir_save):
        os.makedirs(dir_save)

    idx_year = 0
    while year < maxyear:

        rest_years = int(maxyear) - int(year) +1
        if rest_years <10:
            dy = rest_years

        #start_time = datetime.now()

        # calculate number of days within 10-year period
        # predefine array that gets filled with values
        if str(cal) == '360_day':
            ndays = 360*dy
            var_array = np.zeros((ndays,len(lat1d),len(lon1d)))*np.nan
        elif str(cal) == '365_day':
            ndays = 365*dy
            var_array = np.zeros((ndays,len(lat1d),len(lon1d)))*np.nan
        else:
            d0 = date(year, 1, 1)
            d1 = date(year+dy, 1, 1)
            delta = d1-d0
            ndays = delta.days
            var_array = np.zeros((ndays,len(lat1d),len(lon1d)))*np.nan

        print (year)

        for y in range(0,len(lat1d)):
            for x in range(0,len(lon1d)):
                fname = sorted(glob.glob(dir_ts+param+'_y'+str(y)+'_x'+str(x)+'.nc'))
                # check if file exists
                if fname:
                    ds_tmp = xr.open_dataset(fname[0])[param][idx_year:idx_year+ndays]
                    var_array[:,y,x] = ds_tmp.data
                    #ds_tmp = xr.open_dataset(fname[0]).sel(time=slice(str(year), str(year+9)))
                    #var_array[:,y,x] = ds_tmp[param].data
                    ds_tmp.close()
                else:
                    continue
        
        # check for unrealistic values
        if (param == 'pr') | (param == 'rsds') | (param == 'sfcWind') | (param == 'hurs'):
            var_array[(var_array) < 0.] = 0. 
        if (param == 'hurs'):
            var_array[var_array > 100.] = 100.
        
        var_array[np.isnan(var_array)] = -9999
        #%%      
        # create netcdf file
        savename = dir_save+param+'_'+model_name+'_'+str(year)+'-'+str(year+dy-1)+'.nc'
        print ('saving file as netcdf:'+savename)  
        # check if file already exists - if yes, delete it before writing new file        
        if os.path.exists(savename):
            os.remove(savename)

        # open new netCDF file in write mode:
        dataset = Dataset(savename,'w',format='NETCDF4_CLASSIC')

        # create dimensions
        dataset.createDimension("time",None)
        lat = dataset.createDimension("lat", len(lat1d))
        lon = dataset.createDimension("lon",len(lon1d))
        #bnds = dataset.createDimension("bnds",2)

        # create variables
        times = dataset.createVariable("time","f8",("time",), zlib=True)
        lats = dataset.createVariable("lat","f4",("lat",), zlib=True, fill_value = -9999)
        lons = dataset.createVariable("lon","f4",("lon",), zlib=True, fill_value = -9999)
        #time_bnds = dataset.createVariable("time_bnds", "f4",("time","bnds"))
        var = dataset.createVariable(param, "f4", ("time","lat","lon",), zlib=True, fill_value = -9999)
        crs = dataset.createVariable('crs', 'i', ())

        # add attributes
        times.units = 'days since 1950-01-01T00:00:00Z'
        # get calendar from original file
        times.calendar = cal
        #times.bounds = 'time_bnds'
        times.standard_name = 'time'
        times.long_name = 'time'
        times.axis = 'T'

        lats.units = 'degrees_north'
        lats.long_name = 'latitude'
        lats.standard_name = 'latitude'
        lats._CoordinateAxisType = 'Lat'

        lons.units = 'degrees_east'
        lons.long_name = 'longitude'
        lons.standard_name = 'longitude'
        lons._CoordinateAxisType = 'Lon'
            
        if param == 'pr':
            var.units = 'mm'
            var.long_name = 'total daily precipitation'
            var.stnadard_nbame = 'precipitation amount'
            dataset.title='Total daily precipitation'
            dataset.comment='Total daily precipitation bias corrected (scaled distribution mapping) data of the CORDEX model '+model_name+'. The reference period is 1981-2010, the years 2006-2010 are taken from the corresponding rcp4.5 scenario.'
            #var[:] = var_array
            
        elif param == 'tasmax':
            var.units = 'degree_Celsius'
            var.long_name = 'daily maximum near-surface air temperature'
            var.standard_name = 'air_temperature'
            dataset.title='Daily maximum near-surface air temperature'
            dataset.comment='Daily maximum near-surface air temperature bias corrected (scaled distribution mapping) data of the CORDEX model '+model_name+'. The reference period is 1981-2010, the years 2006-2010 are taken from the corresponding rcp4.5 scenario.'
            #var[:] = var_array + 273.15
            
        elif param == 'tasmin':
            var.units = 'degree_Celsius'
            var.long_name = 'daily minimum near-surface air temperature'
            var.standard_name = 'air_temperature'       
            dataset.title='Daily minimum near-surface air temperature'
            dataset.comment='Daily minimum near-surface air temperature bias corrected (scaled distribution mapping) data of the CORDEX model '+model_name+'. The reference period is 1981-2010, the years 2006-2010 are taken from the corresponding rcp4.5 scenario.'
            #var[:] = var_array + 273.15
            
        elif param == 'rsds':
            var.units = 'W m-2'
            var.long_name = 'surface downwelling shortwave flux'
            var.standard_name = 'surface_downwelling_shortwave_flux_in_air'       
            dataset.title='Daily global radiation'
            dataset.comment='Daily global radiation bias corrected (scaled distribution mapping) data of the CORDEX model '+model_name+'. The reference period is 1981-2010, the years 2006-2010 are taken from the corresponding rcp4.5 scenario.'
            #var[:] = var_array + 273.15
            
        elif param == 'sfcWind':
            var.units = 'm s-1'
            var.long_name = 'daily mean 10-m wind speed'
            var.standard_name = 'wind_speed'       
            dataset.title='Daily mean 10-m wind speed'
            dataset.comment='Daily mean 10-m wind speed bias corrected (scaled distribution mapping) data of the CORDEX model '+model_name+'. The reference period is 1981-2010, the years 2006-2010 are taken from the corresponding rcp4.5 scenario.'
    
        elif param == 'hurs':
            var.units = 'percent'
            var.long_name = 'daily mean relative humidity'
            var.standard_name = 'relative_humidity'       
            dataset.title='Daily mean relative humidity'
            dataset.comment='Daily mean relative humidity bias corrected (scaled distribution mapping) data of the CORDEX model '+model_name+'. The reference period is 1981-2010, the years 2006-2010 are taken from the corresponding rcp4.5 scenario.'
            
        var.grid_mapping = 'latitude_longitude'
            
                
        crs.grid_mapping_name = "latitude_longitude"
        crs.longitude_of_prime_meridian = 0.0 
        crs.semi_major_axis = 6378137.0
        crs.inverse_flattening = 298.257223563
        crs.comment = 'Latitude and longitude on the WGS 1984 datum'

        # write data to netCDF variable
        var[:] = var_array
        lats[:] = lat1d
        lons[:] = lon1d

        # fill in times
        #dates = [datetime(year,1,1)+k*timedelta(days=1) for k in range(var_array.shape[0])]
        if str(cal) == '365_day':
            d = np.arange((year-1950)*365, (year-1950)*365 + ndays)
            dates = num2date(d, 'days since 1950-01-01T00:00:00Z', calendar=cal)                            
        elif str(cal) == '360_day':
            d = np.arange((year-1950)*360, (year-1950)*360 + ndays)
            dates = num2date(d, 'days since 1950-01-01T00:00:00Z', calendar=cal)              
        else:
            dates = [datetime(year,1,1)+k*timedelta(days=1) for k in range(var_array.shape[0])]

        times[:] = date2num(dates, units=times.units, calendar=times.calendar)

        # global attributes
        dataset.modelname = model_name
        
        dataset.project="Climaproof, funded by the Austrian Development Agency (ADA) and co-funded by the United Nations Environmental Programme (UNEP)"
        dataset.method = 'Bias Correction Method: Scaled Distribution Mapping (Switanek et al., 2017, doi.org/10.5194/hess-21-2649-2017)'
        dataset.source = 'Created with the ICC-OBS Tool (Institute of Meteorology, University of Natural Resources and Life Sciences, Vienna, Austria)'
        dataset.conventions = 'CF-1.6'
        
        # close dataset        
        dataset.close()
        #%%

        year = year+dy   
        idx_year = idx_year+ndays

        #print 'reordering and saving '+str(dy)+' years of data took '+str(datetime.now() - start_time)
        
    return
    
