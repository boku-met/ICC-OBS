'''
##########################################################
#####Copyright (c) 2016, Matthew Switanek ################
########## edited by Maria Wind, 2018 ####################
##########################################################
'''

import numpy as np
import scipy.stats as stats
import warnings

def fit_distribution(sort_obs, distr):
     
    if distr == 'gamma':
        dist_fit = stats.gamma.fit(sort_obs)
        cdf_fit = stats.gamma.cdf(sort_obs,*dist_fit)
        
    elif distr == 'weibull':
        dist_fit = stats.weibull_min.fit(sort_obs, floc=min(sort_obs)-.1)
        cdf_fit = stats.weibull_min.cdf(sort_obs,*dist_fit)
        
    elif distr == 'normal':
        dist_fit = [np.mean(sort_obs),np.std(sort_obs)]
        cdf_fit = stats.norm.cdf(sort_obs,*dist_fit)

    return dist_fit, cdf_fit


# before applying function, check if observations and model data have same units!

# absolute bias correction (temperature, humidity, wind)

def sdm_absolute(obs_calib, obs_calib_ts, mod_calib, mod_calib_ts, mod_bcperiod, mod_bcperiod_ts, start_f, end_f, distr = 'normal'):
    '''
    distr   string describing distribution that should be fitted to the data
            options: normal, weibull, gamma
    '''
    
    #This is an empty array that we will fill with the bias corrected values for the future model.
    #If you treat the future period as the calibration data, then you will bias correct the historical model data.
    bc_vals_bcperiod = np.zeros(mod_bcperiod.shape,dtype=np.float32)
    bc_vals_bcperiod[:] = np.nan


    for month in range(1,13):

        #This is the time vector for the calibration period
        ftime1_obs = np.nonzero((obs_calib_ts.month==month) | (obs_calib_ts.month==month%12+1) | (obs_calib_ts.month==(month-2)%12 +1))
        try:
            ftime1_mod = np.nonzero((mod_calib_ts.month==month) | (mod_calib_ts.month==month%12+1) | (mod_calib_ts.month==(month-2)%12 +1))
        except(AttributeError, ValueError):
            mod_calib_month = mod_calib_ts.dt.month.data
            ftime1_mod = np.nonzero((mod_calib_month==month) | (mod_calib_month==month%12+1) | (mod_calib_month==(month-2)%12 +1))
        #The observed and modeled data for the calibration time for the given month.
        obs_curr_calib = obs_calib[ftime1_obs]
        mod_curr_calib = mod_calib[ftime1_mod]
    
        #The next line screens for real data in case of any NaNs.
        #fnd1 = np.nonzero(obs_curr_calib>-300)[0] 
        fnd1 = np.nonzero(~np.isnan(obs_curr_calib))[0]
        
        #The next lines detrend the observed and model data in the calibration period.  
        #detrend observations
        xx1 = np.ones((len(obs_curr_calib),2))
        xx1[:,0] = np.arange(0,len(obs_curr_calib))
        rr1 = np.linalg.lstsq(xx1[fnd1],obs_curr_calib[fnd1], rcond=-1)[0]
        dtrend_obs_calib = obs_curr_calib[fnd1] - np.sum(rr1*xx1[fnd1],axis=1)
        #detrend modeldata
        xx1_mod = np.ones((len(mod_curr_calib),2))
        rr2 = np.linalg.lstsq(xx1_mod,mod_curr_calib,rcond=-1)[0]
        dtrend_mod_calib = mod_curr_calib  - np.sum(rr2*xx1_mod,axis=1)
    
        #Sort the detrended data.
        sort_obs_calib = np.sort(dtrend_obs_calib)
        sort_mod_calib = np.sort(dtrend_mod_calib)
        #Get the mean of the data.
        mean_obs_calib = np.nanmean(obs_curr_calib)
        mean_mod_calib = np.mean(mod_curr_calib)
    
        #Get parameters of the normal distribution and corresponding cdf values
        normfit1, cdf_obs_calib = fit_distribution(sort_obs_calib, distr)
        #normfit1 = [np.mean(sort_obs_calib),np.std(sort_obs_calib)]
        #cdf_obs_calib = stats.norm.cdf(sort_obs_calib,*normfit1)
        normfit2 , cdf_mod_calib = fit_distribution(sort_mod_calib, distr) #[np.mean(sort_mod_calib),np.std(sort_mod_calib)]
        #cdf_mod_calib = stats.norm.cdf(sort_mod_calib,*normfit2)
    
        #Limit the extremity of values.
        ffy_obs_tooextreme = cdf_obs_calib > .9999
        cdf_obs_calib[ffy_obs_tooextreme] = .9999
        ffy_obs_tooextreme = cdf_obs_calib < .0001
        cdf_obs_calib[ffy_obs_tooextreme] = .0001
    
        ffy_mod_tooextreme = cdf_mod_calib > .9999
        cdf_mod_calib[ffy_mod_tooextreme] = .9999
        ffy_mod_tooextreme = cdf_mod_calib < .0001
        cdf_mod_calib[ffy_mod_tooextreme] = .0001
        
        year = start_f
        while year < end_f:      
            
            #This is the time vector for the future period    
            #Time vector for the period of interest (30-year window to correct the middle 10 years)
            try:
                ftime2 = np.nonzero(((mod_bcperiod_ts.month==month) | (mod_bcperiod_ts.month==month%12+1) | (mod_bcperiod_ts.month==(month-2)%12 +1)) & ((mod_bcperiod_ts.year >= year) & (mod_bcperiod_ts.year <= year+29)))                 
            except(AttributeError):
                mod_bc_month = mod_bcperiod_ts.dt.month.data
                mod_bc_year = mod_bcperiod_ts.dt.year.data
                ftime2 = np.nonzero(((mod_bc_month==month) | (mod_bc_month==month%12+1) | (mod_bc_month==(month-2)%12 +1)) & ((mod_bc_year >= year) & (mod_bc_year <= year+29)))                 

            mod_curr_bcsubperiod = mod_bcperiod[ftime2]
            mod_curr_bcsubperiod_ts = mod_bcperiod_ts[ftime2]
            
            #Time vector for period that is bias corrected (middle 10 years of 30 year window)
            if year == start_f:
                start_bc_period = year
                end_bc_period = year+19
            elif year+30 >= end_f:
                start_bc_period = year+10
                end_bc_period   = year+29
            else:
                start_bc_period = year+10
                end_bc_period   = year+19
                
            try:
                ftime3 = np.nonzero((mod_bcperiod_ts.month == month) & (mod_bcperiod_ts.year>=start_bc_period) & (mod_bcperiod_ts.year <= end_bc_period))
                # Time vector for bias corrected period within period of interest
                ftime4 = np.nonzero((mod_curr_bcsubperiod_ts.month == month) & (mod_curr_bcsubperiod_ts.year>=start_bc_period) & (mod_curr_bcsubperiod_ts.year <= end_bc_period))
            except(AttributeError):
                ftime3 = np.nonzero((mod_bc_month == month) & (mod_bc_year >= start_bc_period) & (mod_bc_year <= end_bc_period))
                # Time vector for bias corrected period within period of interest
                mod_curr_year = mod_curr_bcsubperiod_ts.dt.year.data
                ftime4 = np.nonzero((mod_curr_bcsubperiod_ts.dt.month.data == month) & (mod_curr_year >=start_bc_period) & (mod_curr_year <= end_bc_period))
            
            #Detrend the future data and get the mean.
            xx2 = np.ones((len(mod_curr_bcsubperiod),2))
            xx2[:,0] = np.arange(0,len(mod_curr_bcsubperiod))
            mean_mod_bcsubperiod = np.mean(mod_curr_bcsubperiod)
            rr3 = np.linalg.lstsq(xx2,mod_curr_bcsubperiod,rcond=-1)[0]
            dtrend_mod_bcsubperiod = mod_curr_bcsubperiod  - np.sum(rr3*xx2,axis=1)       
        
            #Sort future data and get indices of sorted data.
            sort_mod_bcsubperiod = np.sort(dtrend_mod_bcsubperiod)   
            argsort_mod_bcsubperiod = np.argsort(dtrend_mod_bcsubperiod)
        
            #Fit of future data and corresponding cdf values.
            normfit3, cdf_mod_bcsubperiod = fit_distribution(sort_mod_bcsubperiod, distr)
            #normfit3 = [np.mean(sort_mod_bcsubperiod),np.std(sort_mod_bcsubperiod)]
            #cdf_mod_bcsubperiod = stats.norm.cdf(sort_mod_bcsubperiod,*normfit3)
        
            #Limit the extremity of values.
            ffy_mod_tooextreme =  cdf_mod_bcsubperiod > .9999
            cdf_mod_bcsubperiod[ffy_mod_tooextreme] = .9999
            ffy_mod_tooextreme =  cdf_mod_bcsubperiod < .0001
            cdf_mod_bcsubperiod[ffy_mod_tooextreme] = .0001
        
            #Scaling between the future and past distributions.
            if distr == 'gamma':
                distribution_scaling = (stats.gamma.ppf(cdf_mod_bcsubperiod,*normfit3) - stats.gamma.ppf(cdf_mod_bcsubperiod,*normfit2)) \
                    * (normfit1[1]/normfit2[1]) 
            elif distr == 'weibull':
                distribution_scaling = (stats.weibull_min.ppf(cdf_mod_bcsubperiod,*normfit3) - stats.weibull_min.ppf(cdf_mod_bcsubperiod,*normfit2)) \
                    * (normfit1[1]/normfit2[1]) 
            elif distr == 'normal':
                distribution_scaling = (stats.norm.ppf(cdf_mod_bcsubperiod,*normfit3) - stats.norm.ppf(cdf_mod_bcsubperiod,*normfit2)) \
                    * (normfit1[1]/normfit2[1]) 
        
            #Check to see if lengths of future model and observed period are different. 
            if len(cdf_obs_calib) != len(cdf_mod_bcsubperiod):
                xx1 = np.linspace(1, len(cdf_obs_calib),len(cdf_obs_calib))
                xx1new = np.linspace(1, len(cdf_obs_calib),len(cdf_mod_bcsubperiod))
                cdf_obs_calib2 = np.interp(xx1new, xx1, cdf_obs_calib)                                
            else:
                cdf_obs_calib2 = cdf_obs_calib[:]
               
            #Check to see if lengths of future and past model period are different.  
            if len(cdf_mod_calib) != len(cdf_mod_bcsubperiod):
                xx1 = np.linspace(1, len(cdf_mod_calib),len(cdf_mod_calib))
                xx1new = np.linspace(1, len(cdf_mod_calib),len(cdf_mod_bcsubperiod))
                cdf_mod_calib2 = np.interp(xx1new, xx1, cdf_mod_calib)                                
            else:
                cdf_mod_calib2 = cdf_mod_calib[:]
              
            #Get the scaling of the extremity of values. 
            if distr == 'normal':
                #split tails to reflect two-tailed nature of a normal distribution
                cdf_obs_calib_split_tails = cdf_obs_calib2-.5
                cdf_obs_calib_split_tails_signs = np.sign(cdf_obs_calib_split_tails)
                cdf_mod_calib_split_tails = cdf_mod_calib2-.5
                cdf_mod_bcsubperiod_split_tails = cdf_mod_bcsubperiod-.5     
                cdf_obs_calib_period = 1./(.5-np.abs(cdf_obs_calib_split_tails))
                cdf_mod_calib_period = 1./(.5-np.abs(cdf_mod_calib_split_tails))
                cdf_mod_bcsubperiod_period = 1./(.5-np.abs(cdf_mod_bcsubperiod_split_tails))
                cdf_ratio_bcsubperiod2calib = cdf_mod_bcsubperiod_period / cdf_mod_calib_period 
                cdf2use = np.sort(cdf_obs_calib_split_tails_signs*np.abs(.5-1./np.maximum(1,cdf_obs_calib_period * cdf_ratio_bcsubperiod2calib))+.5) 
            else:
                cdf_obs_calib_period = 1./(1-cdf_obs_calib2)
                cdf_mod_calib_period = 1./(1-cdf_mod_calib2)
                cdf_mod_bcsubperiod_period = 1./(1-cdf_mod_bcsubperiod)
                cdf2use = np.sort(1-1./np.maximum(1,cdf_obs_calib_period*(cdf_mod_bcsubperiod_period/cdf_mod_calib_period)))
            
            #Need realistic values.
            ffnonreal = cdf2use > .9999
            cdf2use[ffnonreal] = .9999
            ffnonreal = cdf2use < .0001
            cdf2use[ffnonreal] = .0001
            
            #Get initial bias corrected values.
            if distr == 'gamma':
                bcvals_bcsubperiod = stats.gamma.ppf(cdf2use,*normfit1) + distribution_scaling
            elif distr == 'weibull':
                bcvals_bcsubperiod = stats.weibull_min.ppf(cdf2use,*normfit1) + distribution_scaling
            elif distr == 'normal':
                bcvals_bcsubperiod = stats.norm.ppf(cdf2use,*normfit1) + distribution_scaling
        
            #Fill into appropriate indices and add the appropriate trend back in.
            unsort_bcvals_bcsubperiod1 = np.zeros((len(mod_curr_bcsubperiod)))
            unsort_bcvals_bcsubperiod1[argsort_mod_bcsubperiod] = bcvals_bcsubperiod[:] - np.mean(bcvals_bcsubperiod) + mean_obs_calib + (mean_mod_bcsubperiod-mean_mod_calib)
            unsort_bcvals_bcsubperiod2 = unsort_bcvals_bcsubperiod1 + np.sum(rr3*xx2,axis=1) - mean_mod_bcsubperiod       
        
            #Fill into appropriate month.
            #Only values of current month and 10-year-period within moving window
            bc_vals_bcperiod[ftime3] = unsort_bcvals_bcsubperiod2[ftime4]
            
            year = year+10
            
        #end while
    #end for
    return bc_vals_bcperiod            

#%% relative bias correction (precipitation, radiation)

def sdm_relative(lower_threshold, obs_calib, obs_calib_ts, mod_calib, mod_calib_ts, mod_bcperiod, mod_bcperiod_ts, start_f, end_f, distr = 'gamma'):

    warnings.filterwarnings("ignore")
    obs_calib[obs_calib < lower_threshold] = 0
    mod_calib[mod_calib < lower_threshold] = 0
    mod_bcperiod[mod_bcperiod < lower_threshold] = 0
    warnings.resetwarnings()
    #We will bias correct by a moving seasonal window (use Jan, Feb & Mar to correct Feb ...)

    #This is an empty array that we will fill with the bias corrected values for the future model.
    #If you treat the future period as the calibration data, then you will bias correct the historical model data.
    bc_vals_bcperiod = np.zeros(mod_bcperiod.shape,dtype=np.float32)
    bc_vals_bcperiod[:] = np.nan

    for month in range(1,13): 
        
        #This is the time vector for the calibration period
        #ftime1 = np.nonzero(obs_calib_ts.month==month)
        ftime1_obs = np.nonzero((obs_calib_ts.month==month) | (obs_calib_ts.month==month%12+1) | (obs_calib_ts.month==(month-2)%12 +1))
        try:
            ftime1_mod = np.nonzero((mod_calib_ts.month==month) | (mod_calib_ts.month==month%12+1) | (mod_calib_ts.month==(month-2)%12 +1))
        except(AttributeError, ValueError):
            mod_calib_month = mod_calib_ts.dt.month.data
            ftime1_mod = np.nonzero((mod_calib_month==month) | (mod_calib_month==month%12+1) | (mod_calib_month==(month-2)%12 +1))
    
        #The observed and modeled data for the calibration time for the given period.
        obs_curr_calib = obs_calib[ftime1_obs]
        mod_curr_calib = mod_calib[ftime1_mod]
    
        #Sort data.  
        sort_obs_calib = np.sort(obs_curr_calib)
        sort_mod_calib = np.sort(mod_curr_calib)
        #Find where the data is > 0.
        fnd1 = np.nonzero(sort_obs_calib > 0)[0]
        fnd2 = np.nonzero(sort_mod_calib > 0)[0]  
        #Fit gamma distributions or more appropriate distributions (e.g., Weibull)
        #and find corresponding cdf values with the fitted distributions.
        gamfit1, cdf_obs_calib = fit_distribution(sort_obs_calib[fnd1], distr)
        gamfit2, cdf_mod_calib = fit_distribution(sort_mod_calib[fnd2], distr)
#        gamfit1 = stats.gamma.fit(sort_obs_calib[fnd1],floc=lower_threshold-.1)
#        gamfit2 = stats.gamma.fit(sort_mod_calib[fnd2],floc=lower_threshold-.1)
#        cdf_obs_calib = stats.gamma.cdf(sort_obs_calib[fnd1],*gamfit1)
#        cdf_mod_calib = stats.gamma.cdf(sort_mod_calib[fnd2],*gamfit2)        
        
        #Limit the extremity of the values.
        ffinf = cdf_obs_calib > .999999
        cdf_obs_calib[ffinf] = .999999
        ffinf = cdf_mod_calib > .999999
        cdf_mod_calib[ffinf] = .999999
    
        #The "future" model data.
        year = start_f
        while year < end_f:              
            
            #This is the time vector for the future period    
            #Time vector for the period of interest (30-year window to correct the middle 10 years)
            try:
                ftime2 = np.nonzero(((mod_bcperiod_ts.month==month) | (mod_bcperiod_ts.month==month%12+1) | (mod_bcperiod_ts.month==(month-2)%12 +1)) & ((mod_bcperiod_ts.year >= year) & (mod_bcperiod_ts.year <= year+29)))                 
            except(AttributeError):
                mod_bc_month = mod_bcperiod_ts.dt.month.data
                mod_bc_year = mod_bcperiod_ts.dt.year.data
                ftime2 = np.nonzero(((mod_bc_month==month) | (mod_bc_month==month%12+1) | (mod_bc_month==(month-2)%12 +1)) & ((mod_bc_year >= year) & (mod_bc_year <= year+29)))                 
               
            mod_curr_bcsubperiod = mod_bcperiod[ftime2]
            mod_curr_bcsubperiod_ts = mod_bcperiod_ts[ftime2]
            
            #Time vector for period that is bias corrected (middle 10 years of 30 year window)
            if year == start_f:
                start_bc_period = year
                end_bc_period = year+19
            elif year+30 >= end_f:
                start_bc_period = year+10
                end_bc_period   = year+29
            else:
                start_bc_period = year+10
                end_bc_period   = year+19

            try:
                ftime3 = np.nonzero((mod_bcperiod_ts.month == month) & (mod_bcperiod_ts.year>=start_bc_period) & (mod_bcperiod_ts.year <= end_bc_period))
                # Time vector for bias corrected period within period of interest
                ftime4 = np.nonzero((mod_curr_bcsubperiod_ts.month == month) & (mod_curr_bcsubperiod_ts.year>=start_bc_period) & (mod_curr_bcsubperiod_ts.year <= end_bc_period))
            except(AttributeError):
                ftime3 = np.nonzero((mod_bc_month == month) & (mod_bc_year >= start_bc_period) & (mod_bc_year <= end_bc_period))
                # Time vector for bias corrected period within period of interest
                mod_curr_year = mod_curr_bcsubperiod_ts.dt.year.data
                ftime4 = np.nonzero((mod_curr_bcsubperiod_ts.dt.month.data == month) & (mod_curr_year >=start_bc_period) & (mod_curr_year <= end_bc_period))
           
            #Sort data. #Find where the data is > 0. 
            sort_mod_bcsubperiod = np.sort(mod_curr_bcsubperiod)
            fnd3 = np.nonzero(sort_mod_bcsubperiod > 0)[0]    
            #Find indices of sorted data.                
            fnd_bcsubperiod_argsort = np.argsort(mod_curr_bcsubperiod)
            
            ### If there are too few data points to fit any reasonable distribution, then it does not make sense to bias correct the data. You can adjust these for more robustness.
            if len(fnd1)>=5 and len(fnd2)>=5 and len(fnd3)>=5:
            
                #Fit gamma distribution or more appropriate distribution (e.g., Weibull).
                gamfit3, cdf_mod_bcsubperiod = fit_distribution(sort_mod_bcsubperiod[fnd3], distr)
#                gamfit3 = stats.gamma.fit(sort_mod_bcsubperiod[fnd3],floc=lower_threshold-.1)
#                cdf_mod_bcsubperiod = stats.gamma.cdf(sort_mod_bcsubperiod[fnd3],*gamfit3)
    
                #Scaling the frequency of the number of days > threshold. 
                frequency_scaling = (len(fnd3)*1./len(mod_curr_bcsubperiod)) * ((len(fnd1)*1./len(obs_curr_calib)) / (len(fnd2)*1./len(mod_curr_calib)))
                #The number of  days > threshold we expect after bias correction for the future model.         
                bc_number_of_raindays = int(np.round(frequency_scaling * len(mod_curr_bcsubperiod),0))
                #Limit the extremity of the values.
                ffinf = cdf_mod_bcsubperiod > .999999
                cdf_mod_bcsubperiod[ffinf] = .999999
    
                #The scaling between the future and past model distributions.
                if distr == 'gamma':
                    distribution_scaling = stats.gamma.ppf(cdf_mod_bcsubperiod,*gamfit3) / stats.gamma.ppf(cdf_mod_bcsubperiod,*gamfit2)
                elif distr == 'weibull':
                    distribution_scaling = stats.weibull_min.ppf(cdf_mod_bcsubperiod,*gamfit3) / stats.weibull_min.ppf(cdf_mod_bcsubperiod,*gamfit2)
                elif distr == 'normal':
                    distribution_scaling = stats.norm.ppf(cdf_mod_bcsubperiod,*gamfit3) / stats.norm.ppf(cdf_mod_bcsubperiod,*gamfit2)
            
                #Limit the model projected changes to a factor of 10 in either direction.
                fnd_unrealistic_changes = distribution_scaling > 10.
                distribution_scaling[fnd_unrealistic_changes] = 10.
                fnd_unrealistic_changes = distribution_scaling < 0.1
                distribution_scaling[fnd_unrealistic_changes] = 0.1
            
                #Interpolate obs calib cdf to number of future model  days > threshold.
                xx1 = np.linspace(1, len(fnd1),len(fnd1))
                xx1new = np.linspace(1, len(fnd1),len(fnd3))
                cdf_obs_calib2 = np.interp(xx1new, xx1, cdf_obs_calib)
            
                #Interpolate mod calib cdf to number of future model  days > threshold.
                xx1 = np.linspace(1, len(fnd2),len(fnd2))
                xx1new = np.linspace(1, len(fnd2),len(fnd3))
                cdf_mod_calib2 = np.interp(xx1new, xx1, cdf_mod_calib)
    
                #Get the scaling of the extremity of values. 
                if distr == 'normal':
                    #split tails to reflect two-tailed nature of a normal distribution
                    cdf_obs_calib_split_tails = cdf_obs_calib2-.5
                    cdf_obs_calib_split_tails_signs = np.sign(cdf_obs_calib_split_tails)
                    cdf_mod_calib_split_tails = cdf_mod_calib2-.5
                    cdf_mod_bcsubperiod_split_tails = cdf_mod_bcsubperiod-.5     
                    cdf_obs_calib_period = 1./(.5-np.abs(cdf_obs_calib_split_tails))
                    cdf_mod_calib_period = 1./(.5-np.abs(cdf_mod_calib_split_tails))
                    cdf_mod_bcsubperiod_period = 1./(.5-np.abs(cdf_mod_bcsubperiod_split_tails))
                    cdf_ratio_bcsubperiod2calib = cdf_mod_bcsubperiod_period / cdf_mod_calib_period 
                    cdf2use = np.sort(cdf_obs_calib_split_tails_signs*np.abs(.5-1./np.maximum(1,cdf_obs_calib_period * cdf_ratio_bcsubperiod2calib))+.5) 
                else: 
                    cdf_obs_calib_period = 1./(1-cdf_obs_calib2)
                    cdf_mod_calib_period = 1./(1-cdf_mod_calib2)
                    cdf_mod_bcsubperiod_period = 1./(1-cdf_mod_bcsubperiod)
                    cdf2use = np.sort(1-1./np.maximum(1,cdf_obs_calib_period*(cdf_mod_bcsubperiod_period/cdf_mod_calib_period)))
                
                #Initial bias corrected data
                if distr == 'gamma':
                    bcvals_bcsubperiod = stats.gamma.ppf(cdf2use,*gamfit1) * distribution_scaling
                elif distr == 'weibull':
                    bcvals_bcsubperiod = stats.weibull_min.ppf(cdf2use,*gamfit1) * distribution_scaling
                elif distr == 'normal':
                    bcvals_bcsubperiod = stats.norm.ppf(cdf2use,*gamfit1) * distribution_scaling
            
                #Fill the data into the appropriate indices. The correct locations in the time series.
                #Adjust frequency of  days > threshold (linear interpolation)
                unsort_bcvals_bcsubperiod = np.zeros((len(mod_curr_bcsubperiod)))
                if len(fnd3) > bc_number_of_raindays:
                    xx1 = np.linspace(1, len(fnd3),len(fnd3))
                    xx1new = np.linspace(1, len(fnd3),bc_number_of_raindays)
                    bcvals_bcsubperiod = np.interp(xx1new, xx1, bcvals_bcsubperiod)
                    unsort_bcvals_bcsubperiod[fnd_bcsubperiod_argsort[-bc_number_of_raindays:]] = bcvals_bcsubperiod[:]
                elif len(fnd3) > 0 and len(fnd3) <= bc_number_of_raindays:
                    unsort_bcvals_bcsubperiod[fnd_bcsubperiod_argsort[-len(fnd3):]] = bcvals_bcsubperiod[-len(fnd3):]
    
                #Fill the data into the appropriate month
                #Take only values of current month and 10-year-period within moving window
                bc_vals_bcperiod[ftime3] = unsort_bcvals_bcsubperiod[ftime4]
                
            year=year+10 
        #end while
    #end for    
    return bc_vals_bcperiod