import os
import sys
import json
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

def format_veg_params(basinMask,lcData,outVeg,scheme='IGBP'):
    
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    if scheme == 'IGBP':
        attriFile = os.path.join(__location__,'veg_type_attributes_igbp.json')
    elif scheme == 'GLCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_glcc.json')
    elif scheme == 'IPCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_ipcc.json')
    else:
        raise NameError('Land cover classification scheme not supported')
    
    band = 1
    
    with open(attriFile) as data_file:    
        attriData = json.load(data_file)
    
    clsAttributes = attriData['classAttributes']

    infiles = [os.path.join(__location__,basinMask),
               os.path.join(__location__,lcData)]
   
    ds = gdal.Open(infiles[0],GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    mask = BandReadAsArray(b1)
    maskRes = ds.GetGeoTransform()[1]
    ds = None
    b1 = None
    
    ds = gdal.Open(infiles[1],GA_ReadOnly)
    b1 = ds.GetRasterBand(band)
    lccls = BandReadAsArray(b1)  
    clsRes = ds.GetGeoTransform()[1]
    ds = None
    b1 = None
    
    ratio = maskRes/clsRes
    
    cells = 0
    
    vegfile = os.path.join(__location__,outVeg)
    
    if os.path.exists(vegfile)==True:
        os.remove(vegfile)
    
    with open(vegfile, 'a') as f:
        
        cnt = 1
        
        for i in range(mask.shape[0]):
            y1 = int(i*ratio)
            y2 = int(y1+ratio)
            for j in range(mask.shape[1]):
                x1 = int(j*ratio)
                x2 = int(x1+ratio)
                
                tmp = lccls[y1:y2,x1:x2]
                                
                if mask[i,j] == 1:
                    if np.any(tmp>len(clsAttributes)-1)==True:
                        negdx = np.where(tmp>len(clsAttributes)-1)
                        try:
                            tmp[negdx] = np.bincount(tmp[np.where(tmp<len(
                                            clsAttributes)-1)].ravel()).argmax()
                        except ValueError:
                            tmp[negdx] = 0
                    uniqcnt = np.unique(tmp)
                    clscnt = np.bincount(tmp.ravel())
                    
                    if type(uniqcnt).__name__ == 'int':
                        Nveg = 1
                    else:
                        Nveg = len(uniqcnt)
                        
                    f.write('{0} {1}\n'.format(cnt,Nveg))
                    
                    if Nveg != 0:
                        for t in range(uniqcnt.size):
                            vegcls = int(uniqcnt[t])
                            Cv = np.float(clscnt[uniqcnt[t]])/np.float(clscnt.sum())
                            attributes = clsAttributes[vegcls]['properties']
                            rdepth1 = str(attributes['rootd1'])
                            rfrac1 = str(attributes['rootfr1'])
                            rdepth2 = str(attributes['rootd2'])
                            rfrac2 = str(attributes['rootfr2'])
                            rdepth3 = str(attributes['rootd3'])
                            rfrac3 = str(attributes['rootfr3'])
                        
                            f.write('\t{0} {1:.4f} {2} {3} {4} {5} {6} {7}\n'.format(vegcls,Cv,rdepth1,rfrac1,rdepth2,rfrac2,rdepth3,rfrac3))           
                
                    cnt+=1
                        
                    cells+=1
                    
                else:
                   cnt+=1
                

                
    return
    
def main():
    format_veg_params(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    return
    
if __name__ == "__main__":
    main()