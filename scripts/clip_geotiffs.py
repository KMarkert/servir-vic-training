import os
import glob
import datetime
import numpy as np
from osgeo import gdal,osr
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

def find_nearest_idx(arr,var):    
    return (np.abs(arr-var)).argmin()
    
spatialRef = osr.SpatialReference()
spatialRef.SetWellKnownGeogCS('WGS_84')

infolder = '~/data/input/forcing/raw/globalprecip/'

# extent [NNN,WWW,SSS,EEE]
bb = [0.15,34.5,-0.75,36.0]

flist = glob.glob(infolder + '*.tif')

outfolder = '~/data/input/forcing/raw/persiann/'

for i in range(len(flist)):
    
    inname = flist[i].split('/')[-1]
    inyr = inname[6:10]
    inmon = inname[10:12]
    inday = inname[12:14]
    
    indate = datetime.date(int(inyr),int(inmon),int(inday))
    
    tt = indate.timetuple()
    
    filename = 'persiann-css.{0}.{1}.NyandoBasin.tif'.format(indate.year,tt.tm_yday)
    
    # read in dataset
    ds = gdal.Open(flist[i])
    
    #get geographic information
    if i == 0:
        gt = ds.GetGeoTransform()
        x_size = ds.RasterXSize
        y_size = ds.RasterYSize
        lon0 = gt[0] 
        lon1 = gt[0] + (x_size*gt[1])
        lat0 = gt[3] + (y_size*gt[-1])
        lat1 = gt[3]
        
        lats = np.arange(lat0,lat1+abs(gt[-1]),abs(gt[-1]))
        lons = np.arange(lon0,lon1+abs(gt[0]),abs(gt[1]))
        
        
    band = ds.GetRasterBand(1)
    dataIn = BandReadAsArray(band)
    
    ds = None
    band = None
    
    ## do the clipping
    img_bb = [find_nearest_idx(lats,bb[0]),find_nearest_idx(lons,bb[1]),
              find_nearest_idx(lats,bb[2]),find_nearest_idx(lons,bb[3])]
    
    dataOut = dataIn[img_bb[2]:img_bb[0],img_bb[1]:img_bb[3]]
    
    #save out dataset
    outpath = outfolder+filename
    print outpath
    
    drv = gdal.GetDriverByName('GTiff')
    
    dsOut = drv.Create(outpath, dataOut.shape[1], dataOut.shape[0], 1, gdal.GDT_Float32)
    dsOut.SetGeoTransform((lons[img_bb[1]], gt[1], 0, lats[img_bb[0]], 0, gt[-1]))
    dsOut.SetProjection(spatialRef.ExportToWkt())
    
    
    band = dsOut.GetRasterBand(1)
    band.WriteArray(dataOut)
    
    dsOut = None
    band  = None
    
