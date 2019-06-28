# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 09:03:27 2018

@author: mariaw
"""

#import statsmodels.api as sm
#import xarray as xr
import numpy as np
#import matplotlib.pyplot as plt
#import glob
from datetime import datetime
import pandas as pd 
import string as s
import glob


def dms2dd(degrees, minutes, seconds):
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60);
    return dd;

#*****************************************************
# read station data
#*****************************************************
def read_stationdata(param, data_dir, fn_metadata, start_y=1981, end_y=2010):
    '''
    Function for nicely formated stationdata (in agreement with the template)
    
    Reads the stationdata and metadata from a .csv file and writes 
    it to a pandas dataframe.
    
    
    INPUT:
    
    param:       name of the parameter (e.g. 'pr')
    
    data_dir:    directory where the stationdata is stored
    
    fn_metadata: filename of the file containing the metadata
    
    start_y:     start year of correction period
    
    end_y:       end year of correction period 
    
    '''   

    data_dir = data_dir+'/'
   
    start_y = str(start_y)
    end_y = str(end_y)
    
    # set missing value
    missval = np.nan    
    
    # read metadata to dataframe
    #station_metadata = np.genfromtxt(data_dir+fn_metadata, delimiter=',', skip_header=1, dtype=[('station_nr','<i8'),('name', 'S19'), (param, 'f8')])    
    try:
        df_metadata = pd.read_csv(data_dir+fn_metadata, na_values={'NaN', 'nan', 'NA', 'na', -99.9, -999, -999.9, -9999})
    except(IOError):
        try:
            df_metadata = pd.read_csv(fn_metadata, na_values={'NaN', 'nan', 'NA', 'na', -99.9, -999, -999.9, -9999})
        except(IOError):
            print('Metadata file not found!')
            return
    
    # create empty dataframe that gets filled with stationdata
    dates = pd.date_range(start_y+'-01-01', end_y+'-12-31', freq = '1D')
    df_st = pd.DataFrame(index = dates, columns = df_metadata.stationnr)#[station_metadata[:,0].astype(int)])
    
    try:
        for stat_nr in df_metadata.stationnr:
            filename = glob.glob(data_dir+str(int(stat_nr))+'_'+param+'.*')
            #data_st = np.genfromtxt(data_dir+str(int(stat_nr))+'_'+param+'.csv', delimiter=',', skip_header=1, dtype=[('date', 'S19'), (param, 'f8')])
#            try:
            data_st = pd.read_csv(filename[-1], index_col=0, na_values={'NaN', 'nan', 'NA', 'na', -99.9, -999, -999.9, -9999})
#            except(IOError, IndexError):
#                print('Exception in stationdata')
#                df_metadata = df_metadata.drop(index=df_metadata.loc[df_metadata.stationnr == stat_nr].index)
#                continue
            #find an fill missing dates with nans        
            idx = pd.date_range(data_st.index[0], data_st.index[-1], freq = '1D')
            data_st = data_st.reindex(idx, fill_value = missval)
            
            # select timeframe 1981-2010        
            st_30yr = data_st[start_y+'-01-01':end_y+'-12-31']
            # fill predefined dataframe with data    
            df_st[stat_nr] = st_30yr   
    except(IOError):
        for stat_name in df_metadata.name:
            print (stat_name)
            stat_nr = df_metadata.loc[df_metadata.name == stat_name].stationnr.values
            #data_st = np.genfromtxt(data_dir+str(int(stat_nr))+'_'+param+'.csv', delimiter=',', skip_header=1, dtype=[('date', 'S19'), (param, 'f8')])
            try:
                filename = data_dir+s.lower(stat_name)+'_'+param+'.*'
                data_st = pd.read_csv(glob.glob(filename)[0], index_col=0, na_values={'NaN', 'nan', 'NA', 'na', -99.9, -999, -999.9, -9999})
            except(IOError, IndexError):
                try: 
                    filename2 =data_dir+param+'_'+s.lower(stat_name)+'.*'
                    data_st = pd.read_csv(glob.glob(filename2)[0], index_col=0, na_values={'NaN', 'nan', 'NA', 'na', -99.9, -999, -999.9, -9999})
                except (IOError, IndexError):
                    print('File '+filename+' or '+filename2+' does not exist!')
                    #remove column in metadata
                    df_metadata = df_metadata.drop(index=df_metadata.loc[df_metadata.name == stat_name].index)
                    continue
           
            #find an fill missing dates with nans        
            idx = pd.date_range(data_st.index[0], data_st.index[-1], freq = '1D')
            data_st = data_st.reindex(idx, fill_value = missval)
            
            # select timeframe 1981-2010        
            st_30yr = data_st[start_y+'-01-01':end_y+'-12-31']
            # fill predefined dataframe with data    
            df_st[stat_nr] = st_30yr   
    
    # return the dataframes    
    return df_st, df_metadata


#list of all station data files 
def read_stationdata_testcase(flist):
    '''
    This function was especially created for the testcase during the Summer School.
    The Stationdata used here is NOT in a very nice format. :(
    '''

    # number of rows to skip (header)
    skiprows = 25
    # define dates that should be included in the final timeseries
    dates = pd.date_range('1981-01-01 07:00:00','2010-12-31 07:00:00', freq='1D')
    
    # predefine an array that will get filled in the loop
    station_metadata = np.zeros((len(flist),4))
    
    # loop over every file in the directory (stored in the object flist)
    for ff in flist:
        
        # get the stationnumber from the filename
        statnr = ff[-10:-4]
        print statnr
        
        # open header as list (reading line by line)
        with open(ff) as myfile:
            header = [next(myfile) for x in xrange(skiprows)]
        
        # get row with lat/lon information
        matching = [j for j in header if "Exportzeitreihe" in j]
        idx_ll = header.index(matching[0])
        lat_lon_row = header[idx_ll-1]
        cnt = 0
        for i,j in enumerate(lat_lon_row):    
            if j == ';':
                if cnt == 0:
                    id1 = i
                    cnt = cnt+1
                else:
                    id2 = i
        
        # get all the lat and lon data (in degrees, minutes and seconds)
        lon_d = lat_lon_row[id1+1:id1+3]
        lon_m = lat_lon_row[id1+4:id1+6]
        lon_s = lat_lon_row[id1+7:id1+9]
        
        lat_d = lat_lon_row[id2+1:id2+3]
        lat_m = lat_lon_row[id2+4:id2+6]
        lat_s = lat_lon_row[id2+7:id2+9]
        
        # convert the lat lon data to decimal degrees
        lon = dms2dd(lon_d, lon_m, lon_s)    
        lat = dms2dd(lat_d, lat_m, lat_s)
        
        # get row with height information
        matching = [k for k in header if "Geographische Koordinaten" in k]
        idx_hgt = header.index(matching[0])
        hgt_row = header[idx_hgt-1]
        hgt = float(hgt_row[hgt_row.index(';')+1:])
        
        # writ data found in the header int an array
        station_metadata[flist.index(ff),:] = [int(statnr), lat, lon, hgt]
            
        # get row where header ends and data starts
        dt_st = [i for i in header if "Werte:" in i]
        data_start = header.index(dt_st[0])+1
        
        # read time column  
        station_dates  = np.genfromtxt(ff, delimiter=';', skip_header=data_start, dtype=[('date', 'S19'), ('value', 'S4')], autostrip=True, usecols=0)
        dates_st = [datetime.strptime(d[0], '%d.%m.%Y %H:%M:%S') for d in station_dates]
        
        # read column with data into array and convert from string to float
        station_data  = np.genfromtxt(ff, delimiter=';', skip_header=data_start, dtype=[('date', 'S19'), ('value', 'S4')], autostrip=True, usecols=1)
        data_st = np.array([s[0].replace(',', '.') for s in station_data])
        # parse strings to numbers, invalid parsing will be set as NaN
        data_st = pd.to_numeric(data_st, errors='coerce')
        
        # write data to pandas series and fill missing values with nan
        st = pd.Series(data_st, index=dates_st)
        idx = pd.date_range(dates_st[0], dates_st[-1], freq='1D')
        st = st.reindex(idx, fill_value=np.nan)
        #dates = [datetime(1971,1,1)+k*timedelta(days=1) for k in range(stationdata.shape[0])]
        
        # select time frame 1981 - 2010
        st_30yr = st['1981-01-01':'2010-12-31']
        
        if ff==flist[0]:
            #create dataframe
            df_st = pd.DataFrame(st_30yr, index = dates, columns = [statnr])
        else:
            df_st[statnr] = st_30yr       
    
    #df_st_mon = df_st.resample('M').sum()
    df_metadata = pd.DataFrame(station_metadata, columns = ['station_nr', 'lat', 'lon', 'height'])   
    
    return(df_st, df_metadata)
        

