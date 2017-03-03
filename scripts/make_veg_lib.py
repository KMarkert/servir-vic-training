import os
import sys
import glob
import json
import warnings
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
from scipy import ndimage

warnings.simplefilter("ignore")

def make_veg_lib(LCFile,LAIFolder,ALBFolder,outVeg,scheme='IGBP'):
    
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    if scheme == 'IGBP':
        attriFile = os.path.join(__location__,'veg_type_attributes_igbp.json')
        waterCls = 0
    elif scheme == 'GLCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_glcc.json')
        waterCls = 12
    elif scheme == 'IPCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_ipcc.json')
        waterCls = 0
    elif scheme == 'RHEAS':
        attriFile = os.path.join(__location__,'veg_type_attributes_rheas.json')
        waterCls = 0
    else:
        raise NameError('Land cover classification scheme not supported')
        
    band = 1
        
    with open(attriFile) as data_file:    
        attriData = json.load(data_file)
    
    clsAttributes = attriData['classAttributes']
    
    ds = gdal.Open(os.path.join(__location__,LCFile),GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    lccls = BandReadAsArray(b1)  
    clsRes = ds.GetGeoTransform()[1]
    ds = None
    b1 = None
        
    veglib = os.path.join(__location__,outVeg)
    
    if os.path.exists(veglib)==True:
        os.remove(veglib)
        
    try:
        GVFFolder
    except:
        GVFFolder=None
    
    laifiles = glob.glob(os.path.join(__location__,LAIFolder,'*.tif'))
    albfiles = glob.glob(os.path.join(__location__,ALBFolder,'*.tif'))
    if GVFFolder != None:
        gvffiles = glob.glob(GVFFolder+'*.tif')
       
    for i in range(12):
        laids = gdal.Open(laifiles[i],GA_ReadOnly)
        b1 = laids.GetRasterBand(band)
        lsRes = laids.GetGeoTransform()[1]
        
        zoomFactor = lsRes / clsRes
        
        laidata = ndimage.zoom(BandReadAsArray(b1),zoomFactor,order=0)
        
        if i == 0:
            laiMon = np.zeros([laidata.shape[0],laidata.shape[1],12])
            albMon = np.zeros([laidata.shape[0],laidata.shape[1],12])
            if GVFFolder != None:
                gvfMon = np.zeros([laidata.shape[0],laidata.shape[1],12])
                
        laiMon[:,:,i] = laidata[:,:]
        
        laids = None
        b1 = None
        
        albds = gdal.Open(albfiles[i],GA_ReadOnly)
        b1 = albds.GetRasterBand(band)
        albdata = ndimage.zoom(BandReadAsArray(b1),zoomFactor,order=0)
        albMon[:,:,i] = albdata[:,:]
        
        albds = None
        b1 = None
        
        if GVFFolder != None:
            gvfds = gdal.Open(gvffiles[i],GA_ReadOnly)
            b1 = gvfds.GetRasterBand(band)
            gvfdata = ndimage.zoom(BandReadAsArray(b1),zoomFactor,order=0)
            gvfMon[:,:,i] = gvfdata[:,:]
            
            gvfds = None
            b1 = None
            
    albMon[np.where(albMon>=1000)] = np.nan
    if GVFFolder != None:
                gvfMon[np.where(gvfMon>=1000)] = np.nan
        
    with open(veglib, 'a') as f:    
        for i in range(len(clsAttributes)):
            
            attributes = clsAttributes[i]['properties']
                    
            if i == waterCls:
                lai = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
                alb = [0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08]
                if GVFFolder != None:
                    gvf = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]         
            else:
                lai = []
                alb = []
                gvf = []
                for j in range(12):
                    clsidx = np.where(lccls==i)
                    
                    laiStep = laiMon[:,:,j]
                    albStep = albMon[:,:,j]
                    
                    lai.append(np.nanmean(laiStep[clsidx])/1000.)
                    alb.append(np.nanmean(albStep[clsidx])/1000.)
                    
                    if GVFFolder != None:
                        gvfStep = gvfMon[:,:,j]
                        gvf.append(np.nanmean(gvfStep[clsidx])/1000.)
            
            overstory = int(attributes['overstory'])
            rarc = str(attributes['rarc'])
            rmin= str(attributes['rmin'])
            vegHeight = float(attributes['h'])
            rgl = str(attributes['rgl'])
            rad_atten = str(attributes['rad_atn'])
            wind_atten = str(attributes['wnd_atn'])
            trunk_ratio = str(attributes['trnk_r'])
        
            rough = 0.123 * vegHeight
            dis = 0.67 * vegHeight
            
            if overstory == 1:
                wind_h = vegHeight+10
            else:
                wind_h = vegHeight+2
            
            comment = str(attributes['classname'])
                
            if GVFFolder != None:
                if i == 0:
                    f.write('#Class\tOvrStry\tRarc\tRmin\tJAN-LAI\tFEB-LAI\tMAR-LAI\tAPR-LAI\tMAY-LAI\tJUN-LAI\tJUL-LAI\tAUG-LAI\tSEP-LAI\tOCT-LAI\tNOV-LAI\tDEC-LAI\tJAN-VEG\tFEB-VEG\tMAR-VEG\tAPR-VEG\tMAY-VEG\tJUN-VEG\tJUL-VEG\tAUG-VEG\tSEP-VEG\tOCT-VEG\tNOV-VEG\tDEC-VEG\tJAN-ALB\tFEB_ALB\tMAR-ALB\tAPR-ALB\tMAY-ALB\tJUN-ALB\tJUL-ALB\tAUG-ALB\tSEP-ALB\tOCT-ALB\tNOV-ALB\tDEC-ALB\tJAN-ROU\tFEB-ROU\tMAR-ROU\tAPR-ROU\tMAY-ROU\tJUN-ROU\tJUL-ROU\tAUG-ROU\tSEP-ROU\tOCT-ROU\tNOV-ROU\tDEC-ROU\tJAN-DIS\tFEB-DIS\tMAR-DIS\tAPR-DIS\tMAY-DIS\tJUN-DIS\tJUL-DIS\tAUG-DIS\tSEP-DIS\tOCT-DIS\tNOV-DIS\tDEC-DIS\tWIND_H\tRGL\trad_atten\twind_atten\ttruck_ratio\tCOMMENT\n')
                    
                f.write('{0}\t{1}\t{2}\t{3}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t{9:.4f}\t{10:.4f}\t{11:.4f}\t{12:.4f}\t{13:.4f}\t{14:.4f}\t{15:.4f}\t{16:.4f}\t{17:.4f}\t{18:.4f}\t{19:.4f}\t{20:.4f}\t{21:.4f}\t{22:.4f}\t{23:.4f}\t{24:.4f}\t{25:.4f}\t{26:.4f}\t{27:.4f}\t{28:.4f}\t{29:.4f}\t{30:.4f}\t{31:.4f}\t{32:.4f}\t{33:.4f}\t{34:.4f}\t{35:.4f}\t{36:.4f}\t{37:.4f}\t{38:.4f}\t{39:.4f}\t{40}\t{40}\t{40}\t{40}\t{40}\t{40}\t{40}\t{40}\t{40}\t{40}\t{40}\t{40}\t{41}\t{41}\t{41}\t{41}\t{41}\t{41}\t{41}\t{41}\t{41}\t{41}\t{41}\t{41}\t{42}\t{43}\t{44}\t{45}\t{46}\t{47}\n'.format(i,
                        overstory,rarc,rmin,lai[0],lai[1],lai[2],lai[3],lai[4],lai[5],lai[6],lai[7],lai[8],lai[9],lai[10],lai[11],gvf[0],gvf[1],gvf[2],gvf[3],gvf[4],gvf[5],gvf[6],gvf[7],gvf[8],gvf[9],gvf[10],gvf[11],alb[0],alb[1],alb[2],alb[3],alb[4],alb[5],alb[6],alb[7],alb[8],alb[9],alb[10],alb[11],rough,dis,wind_h,rgl,rad_atten,wind_atten,trunk_ratio,comment))
                
            else:
                if i == 0:
                    f.write('#Class\tOvrStry\tRarc\tRmin\tJAN-LAI\tFEB-LAI\tMAR-LAI\tAPR-LAI\tMAY-LAI\tJUN-LAI\tJUL-LAI\tAUG-LAI\tSEP-LAI\tOCT-LAI\tNOV-LAI\tDEC-LAI\tJAN-ALB\tFEB_ALB\tMAR-ALB\tAPR-ALB\tMAY-ALB\tJUN-ALB\tJUL-ALB\tAUG-ALB\tSEP-ALB\tOCT-ALB\tNOV-ALB\tDEC-ALB\tJAN-ROU\tFEB-ROU\tMAR-ROU\tAPR-ROU\tMAY-ROU\tJUN-ROU\tJUL-ROU\tAUG-ROU\tSEP-ROU\tOCT-ROU\tNOV-ROU\tDEC-ROU\tJAN-DIS\tFEB-DIS\tMAR-DIS\tAPR-DIS\tMAY-DIS\tJUN-DIS\tJUL-DIS\tAUG-DIS\tSEP-DIS\tOCT-DIS\tNOV-DIS\tDEC-DIS\tWIND_H\tRGL\trad_atten\twind_atten\ttruck_ratio\tCOMMENT\n')
                    
                f.write('{0}\t{1}\t{2}\t{3}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t{9:.4f}\t{10:.4f}\t{11:.4f}\t{12:.4f}\t{13:.4f}\t{14:.4f}\t{15:.4f}\t{16:.4f}\t{17:.4f}\t{18:.4f}\t{19:.4f}\t{20:.4f}\t{21:.4f}\t{22:.4f}\t{23:.4f}\t{24:.4f}\t{25:.4f}\t{26:.4f}\t{27:.4f}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\t{35}\n'.format(i,
                        overstory,rarc,rmin,lai[0],lai[1],lai[2],lai[3],lai[4],lai[5],lai[6],lai[7],lai[8],lai[9],lai[10],lai[11],alb[0],alb[1],alb[2],alb[3],alb[4],alb[5],alb[6],alb[7],alb[8],alb[9],alb[10],alb[11],rough,dis,wind_h,rgl,rad_atten,wind_atten,trunk_ratio,comment))

def main():
    make_veg_lib(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    return

if __name__ == "__main__":
    main()
