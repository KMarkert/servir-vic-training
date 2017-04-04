# -*- coding: utf-8 -*-
#******************************************************************************
# FILE: rout_vic.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script runoff and baseflow outputs from the VIC model formatted
#           as a netCDF and calculates streamflow at an outlet based on the 
#           Lohmann et al. (1996) routing model
# DEPENDENCIES: numpy, pandas, netCDF4, osgeo (gdal)
#******************************************************************************

# import dependencies
import sys
import datetime
import warnings
import numpy as np
import pandas as pd
import netCDF4 as nc
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
from dd_conversions import *

# set system to ignore simple warnings
warnings.simplefilter("ignore")

def save_Q_timeseries(outQ,outFile,startTime,endTime,daily='True'):
    """
    Function saves streamflow time series to CSV file with date information
    
    Arguments:
        outQ -- streamflow time series to be saved
        outFile -- output file path with CSV extension
        startTime -- beginning date of time series as datetime object
        endTime -- ending date of time series  as datetime object
        
    Keywords:
        daily -- Boolean: set true if time series is daily streamflow,
                 false for monthly streamflow. Default: True
        
    Returns:
        n/a
        
    Notes:
        Assumes that streamflow is monthly mean streamflow if not daily
        Does not return a value but writes an output file
        
    """
    
    # get year information
    years = np.arange(startTime.year,endTime.year+1)
    tdelta = endTime-startTime # find offset
    
    # set blank lists for data to be input
    yrs = []
    mon = []
    day = []
    q = []
    
    # check if data is daily
    if daily == 'True':
    
        # loop over each day
        for i in range(tdelta.days+1):
            date = startTime + datetime.timedelta(i) # get date information
            newDate = date.strftime('%Y%m%d') # convert date to readable text
            # start appending data
            yrs.append(newDate[:4])
            mon.append(newDate[4:6])
            day.append(newDate[6:])
            q.append(outQ[i])
    
    # if data is monthly
    else:
        
        cnt = 0 # index counter
        # loop over each year
        for i in range(years.size):
            # if first iteration get the starting month in the year
            if i == 0:
                moff = startTime.month 
                mcnt = 13 - startTime.month
            # else assume a full year's worth of data
            else:
                moff = 1
                mcnt = 12
            # check if ending year and find last month in series
            if years[i] == endTime.year:
                mcnt = endTime.month
            # loop over each month and append data
            for j in range(mcnt):
                yrs.append(years[i])
                mon.append(j+moff)
                day.append(1)
                q.append(outQ[cnt])
                cnt+=1

    # convert data to dictionary
    d = {'Year':yrs,'Month':mon,'Date':day,'Discharge':q}
    
    # make dictionary to data frame
    df = pd.DataFrame(data=d)
    
    # save out dataframe to outfile
    df.to_csv(outFile)
            
    return

def dailyQ_2_monthlyQ(dailyQ,yrs):
    """
    Function converts daily streamflow to monthly mean streamflow
    
    Arguments:
        dailyQ -- daily time series of streamflow from a gauge
        yrs -- number of years that the time series covers 
        
    Keywords:
        n/a
        
    Returns:
        monQ -- monthly mean streamflow for the specified time period
        
    Notes:
        This function assumes that the time series starts on Jan. 1st of a year
        and ends on Dec. 31.
        
    """
    
    # create dictionaries for the indices covering each month for leap and non-leap years
    nidx = {0:['01',0,31],1:['02',31,59],2:['03',59,90],3:['04',90,120],
            4:['05',120,151],5:['06',151,181],6:['07',181,212],7:['08',212,243],
            8:['09',243,273],9:['10',273,304],10:['11',304,334],11:['12',334,364]}
       
    lidx = {0:['01',0,31],1:['02',31,60],2:['03',60,91],3:['04',91,121],
            4:['05',121,152],5:['06',152,182],6:['07',182,213],7:['08',213,244],
            8:['09',244,274],9:['10',274,305],10:['11',305,335],11:['12',335,365]}
    
    # get blank array for monthly streamflow
    monQ = np.zeros(yrs.size*12)

    dt = -365 # counter for days
        
    cnt = 0 # counter for monthly streamflow
    
    # loop over the number of years
    for y in range(yrs.size):
        # get number of days for leap or non-leap year
        if yrs[y]%4 == 0:
            idx = lidx
            dt+=366
        else:
            idx = nidx
            dt+=365
        
        # loop over number of months
        for m in range(12):
            
            midx = idx[mons[m]] # grab monthly indices
            
            # shift monthly indces for the processing year
            t1 = midx[1] + dt 
            t2 = midx[2] + dt
            
            # get mean monthly streamflow
            monQ[cnt] = np.nanmean(dailyQ[t1:t2])
            
            cnt+=1 # plus one to streamflow index counter
            
    return monQ
      

def make_irf(xmask,diff,velo):
    """
    Function for creating the impulse response function (IRF) for a grid cell
    
    Arguments:
        xmask -- longest path for any location within the grid cell to reach
                 the stream channel [m]
        diff -- diffusion value
        velo -- overland flow velocity value
                 
    Returns:
        irf -- normalized impulse response function for a grid cell
    
    Keywords:
        diff -- the diffusivity parameter (default 800) [m^2/s]
        velo -- the flow velocity parameter (default 1.5) [m/s]
        
    """
    step = np.arange(0,48) #number of time steps
    ir = np.zeros(step.size) # blank array for impulse response
    t = 0.0 #time initialization in seconds
    
    for i in step:
        t = t + 3600. #time step of one day (in seconds)
        pot = (velo * t - xmask)**2 / (4.0 * diff * t) # function name???
        if pot > 69:
            h = 0.0
        else:
            #impulse response function
            h = 1.0/(2.0*np.sqrt(np.pi*diff)) * (xmask / (t**1.5)) * np.exp(-pot)
        
        ir[i] = h # pass data to array element
        
    irf = ir / ir.sum() #normalize the output IRF
    
    return irf
    
def make_uh(xmask_val,uh_box,diff,velo):
    """
    Function for creating the unit hydrograph for a grid cell
    
    Arguments:
        xmask_val -- longest path for any location within the grid cell to reach
                     the stream channel [m]
        uh_box -- the monthly hydrograph (values 0-1)
        diff -- diffusion value
        velo -- overland flow velocity value
        
    Keywords:
        none
            
    Returns:
        uh -- UH hydrograph for the grid cell based on the IRF
    
    """    
    irf = make_irf(xmask_val,diff,velo) #get the IRF of the grid cell
    
    #defining constants and setting up 
    le = 48
    ke = 12
    uh_day = 96
    tmax = uh_day * 24
    uh_s = np.zeros(ke+uh_day)
    uh_daily = np.zeros(uh_day)
    fr = np.zeros([tmax,2])
    if uh_box.sum() != 1.:
        #normalize hydrograph if not already
        uh_box = uh_box.astype(np.float) / uh_box.sum()
    
    for i in range(24):
        fr[i,0] = 1. / 24.
    
    for t in range(tmax):
        for l in range(le):
            if t-l > 0:
                fr[t,1] = fr[t,1] + fr[t-l,0] * irf[l]
                
        tt = (t+23) / 24.
        uh_daily[tt-1] = uh_daily[tt-1]+fr[t,1]
        
    # for each time period find the UH based on the IRF
    for k in range(ke):
        for u in range(uh_day):
            uh_s[k+u-1] = uh_s[k+u-1] + uh_box[k] * uh_daily[u]
            
    uh = uh_s / uh_s.sum() #normalize the ouput unit hydrograph
    
    return uh
    
def rout_vic(uhfile,catchfile,roFile,bfFile,outfile,stime,etime,daily=False,
             velocity=1.5,diffusion=800):

    band = 1 #constant value for extracting GIS data
    
    try:
        # read VIC outputs in netCDF format
        # runoff data
        ronc = nc.Dataset(roFile)
        rovar = ronc.variables['runoff']
        # baseflow data
        bfnc = nc.Dataset(bfFile)
        bfvar = bfnc.variables['base']
    except IOError:
        raise IOError('Input flux netCDF files do not exist, check file paths')
    
    # get netCDF time information
    timeunits = bfnc.variables['time'].units.split(' ')[-1].split('-')
    
    # find when netCDF time series start 
    begintime = datetime.date(int(timeunits[0]),int(timeunits[1]),int(timeunits[2]))
    
    #read in stn UH    
    indata = pd.read_csv(uhfile)
                            
    uh_box = np.array(indata.UH)
    
    #read in GIS data for UH generation
    infiles = [catchfile]
            
    try: # to read in raster data
        # loop over each raster and read into memory
        for i in range(len(infiles)):
            #read in GIS raster data into numpy array
            ds = gdal.Open(infiles[i],GA_ReadOnly)
            b1 = ds.GetRasterBand(band)
            var = BandReadAsArray(b1)
            
            if i == 0:
                #create blank 3d array on first pass
                data = np.zeros((var.shape[0],var.shape[1], len(infiles)))
                x_size = ds.RasterXSize
                y_size = ds.RasterYSize
                gt = ds.GetGeoTransform()
                
            data[:,:,i] = var[:,:] # pass data
            
            # Flush
            ds = None
            b1 = None
            
    # if not working, raise error
    except:
        raise IOError('Raster file input error, check that all paths are correct')
    
    # get time information
    t1 = datetime.date(int(stime[:4]),int(stime[4:6]),int(stime[6:]))
    t2 = datetime.date(int(etime[:4]),int(etime[4:6]),int(etime[6:]))
    
    tdelta = t2-t1 # time offset from start to end
    
    t1off = (t1-begintime).days # beginning offset from start of netCDF series
    t2off = (t2-begintime).days+1 # ending offset from start of netCDF series
    
    days = tdelta.days+1 # number of days to run
    yrs = np.arange(int(stime[:4]),int(etime[:4])+1) # number of years to run
                
    idx = np.where(data[:,:,0]>0) #get the cells of the watershed
    
    lfactor = 1E3 #convert mm to m
    afactor = 1E6 #convert km^2 to m^2
    tfactor = 86400. #convert day to second
    
    #x_min, gridSize, 0, y_max, 0, -gridSize
    
    # create arrays of latitude and longitude points
    lons = np.linspace(gt[0],gt[0]+(x_size*gt[1]),x_size)
    lats = np.linspace(gt[3]-(y_size*gt[5]),gt[3],y_size)    
        
    # create blank arrays for ...
    basinUH = np.zeros((idx[0].size,days)) # ...basin wide UH
    gridUH = np.zeros((days+108,108)) # grid based UH
    
    # calculate average resolution for pixels
    res = np.mean([abs(gt[1]),abs(gt[5])])
        
    # iterate over each grid cell contributing to the basin outlet
    for i in range(idx[0].size):
        # convert decimal degrees to meters at lat/lon point
        area_vals = dd2meters([lats[idx[0][i]],lons[idx[1][i]]],res)
        
        # calculate conversion factor for m^3/s
        factor = (data[idx[0][i],idx[1][i],0] * ((area_vals[0]*area_vals[1]))) / (tfactor*lfactor)
        
        # calculate theoretical longest distance from a point in a pixel to stream channel
        xmask_val = ((area_vals[0]*2)+(area_vals[1]*2))/2.
                
        # make UH for grid cell
        uh = make_uh(xmask_val,uh_box,float(diffusion),float(velocity))
        
        # grab runoff and baseflow time series for grid cell
        runoff = (rovar[t1off:t2off,idx[0][i],idx[1][i]]*factor)
        basefl = (bfvar[t1off:t2off,idx[0][i],idx[1][i]]*factor)
        
        
        try: # check that runoff data is available
            tmp = runoff.count()
            continue
        except AttributeError:
            # iterate over each time period in UH
            for u in range(uh.size):
                # calculate the percent contributions from pixel to outlet
                gridUH[u:days+u,u] = (basefl+runoff)*uh[u]
            
            # sum across time period
            basinUH[i,:] = np.nansum(gridUH[:days,:],axis=1)
    
    # sum across basin
    finalQ = np.nansum(basinUH,axis=0)
    
    # close netCDF files
    ronc.close()
    bfnc.close()
    
    # temporal aggregation check
    if daily != 'True':
        outQ = dailyQ_2_monthlyQ(finalQ,yrs)
    else:
        outQ = finalQ
        
    # save out time series
    save_Q_timeseries(outQ,outfile,t1,t2,daily)
    
    return

def main():
    n_args = len(sys.argv)
    
    # Check user inputs
    if n_args != 9:
        print "Wrong user input"
        print "Script performs the routing model for the VIC model"
        print "usage: python rout_vic.py <gauge unit hydrograph> <fraction file> <runoff data> <baseflow data> <output streamflow file> <start time> <end time> <boolean out daily data>"
        print "Exiting system..."
        sys.exit()
    
    else:
        # pass the command line arguments into function
        rout_vic(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
                sys.argv[6],sys.argv[7],sys.argv[8])
    
    return
    
# Execute the main level program if run as standalone    
if __name__ == "__main__":
    main()
