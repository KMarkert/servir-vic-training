# -*- coding: utf-8 -*-
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

warnings.simplefilter("ignore")

def save_Q_timeseries(outQ,outFile,startTime,endTime,daily='True'):
    
    years = np.arange(startTime.year,endTime.year+1)
    tdelta = endTime-startTime
    
    yrs = []
    mon = []
    day = []
    q = []
    
    if daily == 'True':
    
        for i in range(tdelta.days+1):
            date = startTime + datetime.timedelta(i)
            newDate = date.strftime('%Y%m%d')
            yrs.append(newDate[:4])
            mon.append(newDate[4:6])
            day.append(newDate[6:])
            q.append(outQ[i])
    else:
        
        cnt = 0
        for i in range(years.size):
            if i == 0:
                moff = startTime.month 
                mcnt = 13 - startTime.month
            else:
                moff = 1
                mcnt = 12
            if years[i] == endTime.year:
                mcnt = endTime.month
            for j in range(mcnt):
                yrs.append(years[i])
                mon.append(j+moff)
                day.append(1)
                q.append(outQ[cnt])
                cnt+=1

    d = {'Year':yrs,'Month':mon,'Date':day,'Discharge':q}
    
    df = pd.DataFrame(data=d)
    
    df.to_csv(outFile)
            
    return

def dailyQ_2_monthlyQ(dailyQ,yrs):
    mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    nidx = {'Jan':['01',0,31],'Feb':['02',31,59],'Mar':['03',59,90],'Apr':['04',90,120],
        'May':['05',120,151],'Jun':['06',151,181],'Jul':['07',181,212],'Aug':['08',212,243],
        'Sep':['09',243,273],'Oct':['10',273,304],'Nov':['11',304,334],'Dec':['12',334,364]}
       
    lidx = {'Jan':['01',0,31],'Feb':['02',31,60],'Mar':['03',60,91],'Apr':['04',91,121],
        'May':['05',121,152],'Jun':['06',152,182],'Jul':['07',182,213],'Aug':['08',213,244],
        'Sep':['09',244,274],'Oct':['10',274,305],'Nov':['11',305,335],'Dec':['12',335,365]}
    
    monQ = np.zeros(yrs.size*12)
    #
    #if (yrs[0]-1)%4 == 0:
    #    dt = -366
    #else:
    dt = -365
        
    cnt = 0
    
    for y in range(yrs.size):
        if yrs[y]%4 == 0:
            idx = lidx
            dt+=366
        else:
            idx = nidx
            dt+=365
        
        for m in range(len(mons)):
            
            midx = idx[mons[m]]
        
            t1 = midx[1] + dt
            t2 = midx[2] + dt
            
            monQ[cnt] = np.nanmean(dailyQ[t1:t2])
            
            cnt+=1
            
    return monQ      
      

def make_irf(xmask,diff,velo):
    """
    Function for creating the impulse response function (IRF) for a grid cell
    
    Arguments:
        xmask -- longest path for any location within the grid cell to reach
                 the stream channel [m]
    
    Keywords:
        diff -- the diffusivity parameter (default 800) [m^2/s]
        velo -- the flow velocity parameter (default 1.5) [m/s]
        
    """
    step = np.arange(0,48) #number of time steps
    ir = np.zeros(step.size) # blank array for impulse response
    t = 0.0 #time initialization in seconds
    
    for i in step:
        t = t + 3600. #time step of one day (in seconds)
        pot = (velo * t - xmask)**2 / (4.0 * diff * t) #function name???
        if pot > 69:
            h = 0.0
        else:
            #impulse response function
            h = 1.0/(2.0*np.sqrt(np.pi*diff)) * (xmask / (t**1.5)) * np.exp(-pot)
        
        ir[i] = h
        
    irf = ir / ir.sum() #normalize the output IRF
    
    return irf
    
def make_uh(xmask_val,uh_box,diff,velo):
    """
    Function for creating the unit hydrograph for a grid cell
    
    Arguments:
        xmask_val -- longest path for any location within the grid cell to reach
                     the stream channel [m]
        uh_box -- the monthly hydrograph (values 0-1)
        
    Keywords:
        none
    
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
        
    for k in range(ke):
        for u in range(uh_day):
            uh_s[k+u-1] = uh_s[k+u-1] + uh_box[k] * uh_daily[u]
            
    uh = uh_s / uh_s.sum() #normalize the ouput unit hydrograph
    
    return uh
    
def rout_vic(uhfile,catchfile,roFile,bfFile,outfile,stime,etime,daily=False,
             velocity=1.5,diffusion=800):

    band = 1 #constant value for extracting GIS data
    
    ronc = nc.Dataset(roFile)
    rovar = ronc.variables['runoff']
    
    bfnc = nc.Dataset(bfFile)
    bfvar = bfnc.variables['base']
    
    timeunits = bfnc.variables['time'].units.split(' ')[-1].split('-')
    
    begintime = datetime.date(int(timeunits[0]),int(timeunits[1]),int(timeunits[2]))
    
    #read in stn UH    
    indata = pd.read_csv(uhfile)
                            
    uh_box = np.array(indata.UH)
    
    #read in GIS data for UH generation
    infiles = [catchfile]
            
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
            
        data[:,:,i] = var[:,:] #
            
        ds = None
        b1 = None
    
    t1 = datetime.date(int(stime[:4]),int(stime[4:6]),int(stime[6:]))
    t2 = datetime.date(int(etime[:4]),int(etime[4:6]),int(etime[6:]))
    
    tdelta = t2-t1
    
    t1off = (t1-begintime).days
    t2off = (t2-begintime).days+1
    
    days = tdelta.days+1
    yrs = np.arange(int(stime[:4]),int(etime[:4])+1)
                
    idx = np.where(data[:,:,0]>0) #get the cells of the watershed
    
    lfactor = 1E3 #convert mm to m
    afactor = 1E6 #convert km^2 to m^2
    tfactor = 86400. #convert day to second
    
    #x_min, gridSize, 0, y_max, 0, -gridSize
    
    lons = np.linspace(gt[0],gt[0]+(x_size*gt[1]),x_size)
    lats = np.linspace(gt[3]-(y_size*gt[5]),gt[3],y_size)    
        
    basinUH = np.zeros((idx[0].size,days))
    gridUH = np.zeros((days+108,108))
    
    res = np.mean([abs(gt[1]),abs(gt[5])])
        
    for i in range(idx[0].size):
        area_vals = dd2meters([lats[idx[0][i]],lons[idx[1][i]]],res)
        
        #factor = (data[idx[0][i],idx[1][i],0] * (data[idx[0][i],idx[1][i],2]*afactor)) / (tfactor*lfactor)
        factor = (data[idx[0][i],idx[1][i],0] * ((area_vals[0]*area_vals[1]))) / (tfactor*lfactor)
        xmask_val = ((area_vals[0]*2)+(area_vals[1]*2))/4.
                
        uh = make_uh(xmask_val,uh_box,float(diffusion),float(velocity))
        
        runoff = (rovar[t1off:t2off,idx[0][i],idx[1][i]]*factor)
        basefl = (bfvar[t1off:t2off,idx[0][i],idx[1][i]]*factor)
        
        try:
            tmp = runoff.count()
            continue
        except AttributeError:
            for u in range(uh.size):
                gridUH[u:days+u,u] = (basefl+runoff)*uh[u]
                
            basinUH[i,:] = np.nansum(gridUH[:days,:],axis=1)
        
    finalQ = np.nansum(basinUH,axis=0)
    
    ronc.close()
    bfnc.close()
    
    if daily != 'True':
        outQ = dailyQ_2_monthlyQ(finalQ,yrs)
    else:
        outQ = finalQ
    
    save_Q_timeseries(outQ,outfile,t1,t2,daily)
    
    return

def main():
    rout_vic(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
             sys.argv[6],sys.argv[7],sys.argv[8])
    
    return
    
if __name__ == "__main__":
    main()
