import os
import sys
import warnings
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

warnings.simplefilter("ignore")

def format_snow_params(basinMask,elvHiRes,outSnow,interval):
    
    band = 1

    interval = int(interval)
    
    infiles = [basinMask,elvHiRes]
    
    ds = gdal.Open(infiles[0],GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    mask = BandReadAsArray(b1)
    maskRes = ds.GetGeoTransform()[1]
    ds = None
    b1 = None
    
    ds = gdal.Open(infiles[1],GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    elvhires = BandReadAsArray(b1)
    clsRes = ds.GetGeoTransform()[1] 
    ds = None
    b1 = None

    try:
	Precip
    except:
	Precip = None
    
    if Precip != None:
        ds = gdal.Open(infiles[3],GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        pr = BandReadAsArray(b1)
        #prcRes = ds.GetGeoTransform()[1] 
        ds = None
        b1 = None
    
    elvhires = np.ma.masked_where(elvhires<0,elvhires)
    if Precip != None:
        pr = np.ma.masked_where(pr<0,pr)
    
    clsRatio = int(maskRes/clsRes)
    
    #elvXRatio = int(maskRes/elvYRes)
    #elvYRatio = int(maskRes/elvXRes)
    
    snowfile = outSnow
    
    if os.path.exists(snowfile)==True:
        os.remove(snowfile)    

    nbands = []
    
    with open(snowfile, 'a') as f:
        
        cnt = 1
        
        for i in range(mask.shape[0]):
            cy1 = i*clsRatio
            cy2 = cy1+clsRatio
            #ey1 = i*elvYRatio
            #ey2 = ey1+elvYRatio

            for j in range(mask.shape[1]):
        
                cx1 = j*clsRatio
                cx2 = cx1+clsRatio
                #ex1 = j*elvXRatio
                #ex2 = ex1+elvXRatio
                
                tmp = elvhires[cy1:cy2,cx1:cx2]
		
		minelv = np.min(tmp.astype(int)) - (np.min(tmp.astype(int))%interval)
		maxelv = np.max(tmp.astype(int)) + (np.max(tmp.astype(int))%interval)
                                                
		bands = np.arange(minelv, maxelv+interval,interval)

                if mask[i,j] == 1:

		    bcls = np.zeros_like(tmp)
		    bcls[:,:] = -1

		    for b in range(bands.size-1):
			
                    	bcls[np.where((tmp>=bands[b])&(tmp<bands[b+1]))] = b
		    
                    uniqcnt = np.unique(bcls[np.where(tmp>0)])

		    nbands.append(uniqcnt.size)


	maxbands = max(nbands)

	for i in range(mask.shape[0]):
            cy1 = i*clsRatio
            cy2 = cy1+clsRatio
            #ey1 = i*elvYRatio
            #ey2 = ey1+elvYRatio

            for j in range(mask.shape[1]):
        
                cx1 = j*clsRatio
                cx2 = cx1+clsRatio
                #ex1 = j*elvXRatio
                #ex2 = ex1+elvXRatio
                
                tmp = elvhires[cy1:cy2,cx1:cx2]
		
		minelv = np.min(tmp.astype(int)) - (np.min(tmp.astype(int))%interval)
		maxelv = np.max(tmp.astype(int)) + (np.max(tmp.astype(int))%interval)
                                                
		bands = np.arange(minelv, maxelv+interval,interval)

                if mask[i,j] == 1:

		    bcls = np.zeros_like(tmp)
		    bcls[:,:] = -1

		    for b in range(bands.size-1):
			
                    	bcls[np.where((tmp>=bands[b])&(tmp<bands[b+1]))] = b
		    
                    uniqcnt = np.unique(bcls[np.where(tmp>0)])
                    #clscnt = np.bincount(tmp.ravel())                
                    
                    f.write('{0}\t'.format(cnt))
                    
		    lvls = []

                    for c in range(maxbands):
                        try:
                            idx = np.where(bcls==uniqcnt[c])
                            frac = np.float(idx[0].size) / np.float(bcls[np.where(bcls>=0)].size)
                        except IndexError:
                            frac = 0
                        
			f.write('{0:.4f}\t'.format(frac))
                    
                    for c in range(maxbands):
                        try:
                            idx = np.where(bcls==uniqcnt[c])
                            tmpelv = tmp[cy1:cy2,cx1:cx2]
                            muelv = np.nanmean(tmp[idx])
                            
                        except IndexError:
                            muelv = 0
                            
                        f.write('{0:.4f}\t'.format(muelv))

		    for c in range(maxbands):
			try:
                            idx = np.where(bcls==uniqcnt[c])
                            frac = np.float(idx[0].size) / np.float(bcls[np.where(bcls>=0)].size)
                        except IndexError:
                            frac = 0
                        
			f.write('{0:.4f}\t'.format(frac))

		    f.write('\n')


                cnt += 1

    print 'Number of maximum bands: {0}'.format(maxbands)          

    return
                
def main():
    format_snow_params(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    return  

if __name__ == "__main__":
    main()         
   
