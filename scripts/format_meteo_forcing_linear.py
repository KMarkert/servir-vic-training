import os
import sys
import netCDF4
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
import datetime

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
    
    old = 'YYYY.MM.DD'
                
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
        
                    #for ERA/Corrected PERSIANN data
                    ncdfs = [os.path.join('/home/servir-vic/Documents/bc_training/data/Nyando/corrected_persian/YYYY.MM.DD.tif'),
                            os.path.join(__location__,inpath,"nyando_era_met_data.nc")]
                    
                    if yrs == 0:
                    
                        meteonc = netCDF4.Dataset(ncdfs[1])
                        latnc = meteonc.variables['latitude'];latnc=latnc[:]
                        lonnc = meteonc.variables['longitude'];latnc=lonnc[:]
                        time = meteonc.variables['time'][:]
                        
                        tmax = meteonc.variables['mx2t']
                        tmin = meteonc.variables['mn2t']
                        uwind = meteonc.variables['u10']
                        vwind = meteonc.variables['v10']
                        
                        lons,lats = np.meshgrid(lonnc,latnc)
                        
                        idx = find_nearest_idx(lons,lats,x,y)
                        
                        mettime = meteonc.variables
                        
                        
                        cnt1 = 0
			cnt2 = 1
        
                    meteofile = os.path.join(__location__,outpath,'forcing_{0:.4f}_{1:.4f}'.format(y,x))
                    with open(meteofile, 'a') as f:
                        if years[yrs] % 4 == 0:
                            time = 366
                        else:
                            time = 365
                        for t in range(time):			    
			    date = datetime.date(years[yrs],01,01) + datetime.timedelta(t)

			    yr = date.year
			    mn = date.month
			    if len(str(mn)) < 2:
				mn = '0'+str(mn)
			    dy = date.day
			    if len(str(dy)) < 2:
				dy = '0'+str(dy)

			    new = '{0}.{1}.{2}'.format(yr,mn,dy)
                            
                            prds = gdal.Open(ncdfs[0].replace(old,new),GA_ReadOnly)
                            b1 = prds.GetRasterBand(band)
                            pr = BandReadAsArray(b1)
                            
                            if yrs == 0:
                                prX = prds.RasterXSize
                                prY = prds.RasterYSize
                                gt = prds.GetGeoTransform()
                                lonpr0 = gt[0] 
                                lonpr1 = gt[0] + (prX*gt[1])
                                latpr0 = gt[3] + (prY*gt[-1])
                                latpr1 = gt[3]
                                
                                lonpr = np.linspace(lonpr0,lonpr1,prX)
                                latpr = np.linspace(latpr0,latpr1,prY)
                                
                                prlons,prlats = np.meshgrid(lonpr,latpr)
                                pridx = find_nearest_idx(prlons,prlats,x,y)

                                
                            prval = np.mean(pr[pridx[0],pridx[1]])
                            if prval<0:prval=0.0
                            tmaxval = np.max(tmax[cnt1:cnt2+1,idx[0],idx[1]]-273.15)
                            tminval = np.min(tmin[cnt1:cnt2+1,idx[0],idx[1]]-273.15)
    
                            uwndval = uwind[cnt1:cnt2+1,idx[0],idx[1]].mean()
                            vwndval = vwind[cnt1:cnt2+1,idx[0],idx[1]].mean()
                            
                            windval = np.sqrt(uwndval**2 + vwndval**2)
                            
                            f.write('{0} {1} {2} {3}\n'.format(prval,
                                    tmaxval,tminval,windval))     
                                    
                            cnt1+=2
			    cnt2+=2
			    
			    prds = None
			    b1 = None        
                            

                                    
                meteonc.close()
    
    return
    
def main():

    t1 = datetime.datetime.now()
    format_meteo_forcing(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    dt = datetime.datetime.now()-t1
    print 'Processing time: {0}'.format(dt)
    
    return
    
if __name__ == "__main__":
    main()

=======
import os
import sys
import netCDF4
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
import datetime

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
    
    old = 'YYYY/YYYY.MM.DD'
                
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
        
                    #for ERA/Corrected PERSIANN data
                    ncdfs = [os.path.join('/home/servir-vic/Documents/bc_training/data/Nyando/persian/corrected_persian/YYYY/YYYY.MM.DD.tif'),
                            os.path.join(__location__,inpath,"nyando_era_met_data.nc")]
                    
                    if yrs == 0:
                    
                        meteonc = netCDF4.Dataset(ncdfs[1])
                        latnc = meteonc.variables['latitude'];latnc=latnc[:]
                        lonnc = meteonc.variables['longitude'];latnc=lonnc[:]
                        time = meteonc.variables['time'][:]
                        
                        tmax = meteonc.variables['mx2t']
                        tmin = meteonc.variables['mn2t']
                        uwind = meteonc.variables['u10']
                        vwind = meteonc.variables['v10']
                        
                        lons,lats = np.meshgrid(lonnc,latnc)
                        
                        idx = find_nearest_idx(lons,lats,x,y)
                        
                        mettime = meteonc.variables
                        
                        
                        cnt1 = 0
			cnt2 = 1
        
                    meteofile = os.path.join(__location__,outpath,'forcing_{0:.4f}_{1:.4f}'.format(y,x))
                    with open(meteofile, 'a') as f:
                        if years[yrs] % 4 == 0:
                            time = 366
                        else:
                            time = 365
                        for t in range(time):			    
			    date = datetime.date(years[yrs],01,01) + datetime.timedelta(t)

			    yr = date.year
			    mn = date.month
			    if len(str(mn)) < 2:
				mn = '0'+str(mn)
			    dy = date.day
			    if len(str(dy)) < 2:
				dy = '0'+str(dy)

			    new = '{0}/{0}.{1}.{2}'.format(yr,mn,dy)
                            
                            prds = gdal.Open(ncdfs[0].replace(old,new),GA_ReadOnly)
                            b1 = prds.GetRasterBand(band)
                            pr = BandReadAsArray(b1)
                            
                            if yrs == 0:
                                prX = prds.RasterXSize
                                prY = prds.RasterYSize
                                gt = prds.GetGeoTransform()
                                lonpr0 = gt[0] 
                                lonpr1 = gt[0] + (prX*gt[1])
                                latpr0 = gt[3] + (prY*gt[-1])
                                latpr1 = gt[3]
                                
                                lonpr = np.linspace(lonpr0,lonpr1,prX)
                                latpr = np.linspace(latpr0,latpr1,prY)
                                
                                prlons,prlats = np.meshgrid(lonpr,latpr)
                                pridx = find_nearest_idx(prlons,prlats,x,y)

                                
                            prval = np.mean(pr[pridx[0],pridx[1]])
                            if prval<0:prval=0.0
                            tmaxval = np.max(tmax[cnt1:cnt2+1,idx[0],idx[1]]-273.15)
                            tminval = np.min(tmin[cnt1:cnt2+1,idx[0],idx[1]]-273.15)
    
                            uwndval = uwind[cnt1:cnt2+1,idx[0],idx[1]].mean()
                            vwndval = vwind[cnt1:cnt2+1,idx[0],idx[1]].mean()
                            
                            windval = np.sqrt(uwndval**2 + vwndval**2)
                            
                            f.write('{0} {1} {2} {3}\n'.format(prval,
                                    tmaxval,tminval,windval))     
                                    
                            cnt1+=2
			    cnt2+=2
			    
			    prds = None
			    b1 = None        
                            

                                    
                meteonc.close()
    
    return
    
def main():

    t1 = datetime.datetime.now()
    format_meteo_forcing(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    dt = datetime.datetime.now()-t1
    print 'Processing time: {0}'.format(dt)
    
    return
    
if __name__ == "__main__":
    main()
    
>>>>>>> origin/master:scripts/format_meteo_forcing_persiann.py
