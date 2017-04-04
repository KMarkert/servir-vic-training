#------------------------------------------------------------------------------
# FILE: format_meteo_forcing.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script takes meteorological data from ERA-Interim and PERSIANN
#          and writes a time series of precip, tmax, tmin, and wind for each
#          grid cell to run VIC for
# DEPENDENCIES: numpy, netCDF4, osgeo (gdal)
#------------------------------------------------------------------------------

import os
import sys
import netCDF4
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
from datetime import datetime

def find_nearest_idx(xx,yy,xval,yval):    
    xidx = (np.abs(xx-xval)).argmin()
    yidx = (np.abs(yy-yval)).argmin()
    
    ridx = yidx / xx.shape[1]
    cidx = xidx % xx.shape[1]
        
    return [ridx,cidx]
    
def format_meteo_forcing(basin_mask,inpath,outpath,startyr,endyr):
    
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

    band = 1
    
    infiles = [os.path.join(__location__,basin_mask)]
                
    for i in range(len(infiles)):
        
            ds = gdal.Open(infiles[i],GA_ReadOnly)
            b1 = ds.GetRasterBand(band)
            var = BandReadAsArray(b1)
            
            if i == 0:
                data = np.zeros((ds.RasterYSize,ds.RasterXSize, len(infiles)))
                gt = ds.GetGeoTransform()
                lon0 = gt[0] + (gt[1] / 2.)
                lon1 = gt[0] + (data.shape[1]*gt[1])
                lat0 = gt[3] + (data.shape[0]*gt[-1])
                lat1 = gt[3] + (gt[-1] / 2.)
            
            data[:,:,i] = var[:,:]
                
            ds = None
            b1 = None
            
    lons = np.linspace(lon0,lon1,data.shape[1])
    lats = np.linspace(lat0,lat1,data.shape[0])
    
    xx,yy = np.meshgrid(lons,lats)
    
    yy = np.flipud(yy)
    
    years = np.arange(int(startyr),int(endyr)+1)
    
    mask = data[:,:,0].astype(uint8)
    mask = np.ma.masked_where(mask!=1,mask)
    
    for i in range(yy.shape[0]):
        for j in range(yy.shape[1]):
            if mask.mask[i,j] == False:
                x = xx[i,j]
                y = yy[i,j]
            
                for yrs in range(years.size):
        
                    #for ERA/PERSIANN data
                    ncdfs = [os.path.join(__location__,inpath,'persiann/persiann-css.{0}.NyandoBasin.nc'.format(years[yrs])),
                            os.path.join(__location__,inpath,"nyando_era_met_data.nc")]
                            
                    prnc = netCDF4.Dataset(ncdfs[0])
                    pr = prnc.variables['precip']
                    
                    if yrs == 0:
                        
                        latpr = prnc.variables['latitude'];latpr = latpr[:]
                        lonpr = prnc.variables['longitude'];lonpr = lonpr[:]
                    
                        meteonc = netCDF4.Dataset(ncdfs[1])
                        latnc = meteonc.variables['latitude'];latnc=latnc[:]
                        lonnc = meteonc.variables['longitude'];latnc=lonnc[:]
                        time = meteonc.variables['time'][:]
                        
                        tmax = meteonc.variables['mx2t']
                        tmin = meteonc.variables['mn2t']
                        uwind = meteonc.variables['u10']
                        vwind = meteonc.variables['v10']
                        
                        lons,lats = np.meshgrid(lonnc,latnc)
                        prlons,prlats = np.meshgrid(lonpr,latpr)
                        
                        idx = find_nearest_idx(lons,lats,x,y)
                        pridx = find_nearest_idx(prlons,prlats,x,y)
                        
                        cnt1 = 0
			cnt2 = 1
        
                    meteofile = os.path.join(__location__,outpath,'forcing_{0:.4f}_{1:.4f}'.format(y,x))
                    with open(meteofile, 'a') as f:
                        if years[yrs] % 4 == 0:
                            time = 366
                        else:
                            time = 365
                        for t in range(time):
                                
                            prval = np.mean(pr[t,pridx[0],pridx[1]])
                            tmaxval = np.max(tmax[cnt1:cnt2+1,idx[0],idx[1]]-273.15)
                            tminval = np.min(tmin[cnt1:cnt2+1,idx[0],idx[1]]-273.15)
    
                            uwndval = uwind[cnt1:cnt2+1,idx[0],idx[1]].mean()
                            vwndval = vwind[cnt1:cnt2+1,idx[0],idx[1]].mean()
                            
                            windval = np.sqrt(uwndval**2 + vwndval**2)
                            
                            f.write('{0} {1} {2} {3}\n'.format(prval,
                                    tmaxval,tminval,windval))     
                                    
                            cnt1+=2
			    cnt2+=2               
                            

                    prnc.close()                   
                meteonc.close()
    
    return
    
def main():

    # checking user input
    if len(sys.argv) != 6:
        print "Wrong user input"
        print "usage: python format_meteo_forcing.py <template raster> <met data input path> <forcing file outpath> <start year> <end year>"
        #print "DIR INPUTS SHOULD CONTAIN TRAILING /"
        sys.exit()
        
    else:
        if sys.argv[2][-1] != '/':
            print "Input met data dir should contain trailing '/'"
            print "fixing it for you..."
            sys.argv[2] = sys.argv[2] + "/"
            
        if sys.argv[3][-1] != '/':
            print "Output forcing data dir should contain trailing '/'"
            print "fixing it for you..."
            sys.argv[3] = sys.argv[3] + "/"
            
        # pass system arguments to the function
        t1 = datetime.now()
        format_meteo_forcing(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
        dt = datetime.now()-t1
        print 'Processing time: {0}'.format(dt)
    
    return
    
# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()
    
