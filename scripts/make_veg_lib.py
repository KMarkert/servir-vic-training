#******************************************************************************
# FILE: make_veg_lib.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script takes geotiff data for the vegetation library file and 
#          converts it to the needed text format for the VIC model
# DEPENDENCIES: numpy, scipy, osgeo (gdal)
#******************************************************************************

# import dependencies
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

# set system to ignore simple warnings
warnings.simplefilter("ignore")

def make_veg_lib(LCFile,LAIFolder,ALBFolder,outVeg,scheme='IGBP'):
    
    # define script file path for relative path definitions
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    # get land cover classification scheme and create path to lookup table
    if scheme == 'IGBP':
        attriFile = os.path.join(__location__,'veg_type_attributes_igbp.json')
        waterCls = 0
    elif scheme == 'GLCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_glcc.json')
        waterCls = 12
    elif scheme == 'IPCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_ipcc.json')
        waterCls = 0
    else:
        raise SyntaxError('Land cover classification scheme not supported')
        
    band = 1 # constant variable for reading in data
        
    # open/read veg scheme json file
    with open(attriFile) as data_file:    
        attriData = json.load(data_file)
    
    # pass look up information into variable
    clsAttributes = attriData['classAttributes']
    
    try: # try to read in the raster data
    
        ds = gdal.Open(os.path.join(__location__,LCFile),GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        lccls = BandReadAsArray(b1)  
        clsRes = ds.GetGeoTransform()[1]
        ds = None
        b1 = None
        
        # get list of paths to LAI and albedo data
        laifiles = glob.glob(os.path.join(__location__,LAIFolder,'*.tif'))
        albfiles = glob.glob(os.path.join(__location__,ALBFolder,'*.tif'))
        
        # loop over each month in the year
        for i in range(12):
            # read LAI data
            laids = gdal.Open(laifiles[i],GA_ReadOnly)
            b1 = laids.GetRasterBand(band)
            lsRes = laids.GetGeoTransform()[1] # geotransform of land surface data
            
            zoomFactor = lsRes / clsRes # factor for resampling land surface data
            
            # resample LAI land surface data
            laidata = ndimage.zoom(BandReadAsArray(b1),zoomFactor,order=0)
            
            # if first iteration then create blank arrays to pass data to
            if i == 0:
                laiMon = np.zeros([laidata.shape[0],laidata.shape[1],12])
                albMon = np.zeros([laidata.shape[0],laidata.shape[1],12])
                    
            laiMon[:,:,i] = laidata[:,:] # pass lai data in array
            
            # Flush
            laids = None
            b1 = None
            
            # read albedo data
            albds = gdal.Open(albfiles[i],GA_ReadOnly)
            b1 = albds.GetRasterBand(band)
            
            # resmaple albedo land surface data
            albdata = ndimage.zoom(BandReadAsArray(b1),zoomFactor,order=0)
            
            albMon[:,:,i] = albdata[:,:] # pass albedo data in array
            
            # Flush
            albds = None
            b1 = None
            
    # except not working, give error message
    except AttributeError:
        raise IOError('Raster file input error, check that all paths are correct')
        
    # mask nodata values
    albMon[np.where(albMon>=1000)] = np.nan
    
    # get file path to output file
    veglib = os.path.join(__location__,outVeg)
    
    # check if the output parameter file exists, if so delete it
    if os.path.exists(veglib)==True:
        os.remove(veglib)
        
    try: # try to write output veg parameter file
    
        # open output file for writing
        with open(veglib, 'w') as f:
            # loop over each class
            for i in range(len(clsAttributes)):
                
                # get attributes
                attributes = clsAttributes[i]['properties']
                        
                if i == waterCls: # set default values for water
                    lai = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
                    alb = [0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08]
        
                else: # grab lai and albedo data from rasters
                    # set blank lists
                    lai = []
                    alb = []
                    # loop over each month
                    for j in range(12):
                        # find which elements are equal to the class
                        clsidx = np.where(lccls==i)
                        
                        # grab the month time slice
                        laiStep = laiMon[:,:,j]
                        albStep = albMon[:,:,j]
                        
                        # grab the data at location
                        lai.append(np.nanmean(laiStep[clsidx])/1000.)
                        alb.append(np.nanmean(albStep[clsidx])/1000.)
                
                # grab other attributes from lookup table
                overstory = int(attributes['overstory']) # overstory value
                rarc = str(attributes['rarc']) # veg architectural resistance
                rmin= str(attributes['rmin']) # veg minimum stomatal resistance
                vegHeight = float(attributes['h']) # veg height
                rgl = str(attributes['rgl']) # Minimum incoming shortwave radiation for ET
                rad_atten = str(attributes['rad_atn']) # radiation attenuation factor
                wind_atten = str(attributes['wnd_atn']) # wind speed attenuation
                trunk_ratio = str(attributes['trnk_r']) # ratio of total tree height that is trunk
             
                rough = 0.123 * vegHeight # vegetation roughness length 
                dis = 0.67 * vegHeight # vegetation displacement height 
                
                # adjust wind height value if overstory is true
                if overstory == 1:
                    wind_h = vegHeight+10
                else:
                    wind_h = vegHeight+2
            
                comment = str(attributes['classname']) # grab class name
                    
                if i == 0: # write header information
                    f.write('#Class\tOvrStry\tRarc\tRmin\tJAN-LAI\tFEB-LAI\tMAR-LAI\tAPR-LAI\tMAY-LAI\tJUN-LAI\tJUL-LAI\tAUG-LAI\tSEP-LAI\tOCT-LAI\tNOV-LAI\tDEC-LAI\tJAN-ALB\tFEB_ALB\tMAR-ALB\tAPR-ALB\tMAY-ALB\tJUN-ALB\tJUL-ALB\tAUG-ALB\tSEP-ALB\tOCT-ALB\tNOV-ALB\tDEC-ALB\tJAN-ROU\tFEB-ROU\tMAR-ROU\tAPR-ROU\tMAY-ROU\tJUN-ROU\tJUL-ROU\tAUG-ROU\tSEP-ROU\tOCT-ROU\tNOV-ROU\tDEC-ROU\tJAN-DIS\tFEB-DIS\tMAR-DIS\tAPR-DIS\tMAY-DIS\tJUN-DIS\tJUL-DIS\tAUG-DIS\tSEP-DIS\tOCT-DIS\tNOV-DIS\tDEC-DIS\tWIND_H\tRGL\trad_atten\twind_atten\ttruck_ratio\tCOMMENT\n')
                    
                # write the land surface parameterization data
                f.write('{0}\t{1}\t{2}\t{3}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t{9:.4f}\t{10:.4f}\t{11:.4f}\t{12:.4f}\t{13:.4f}\t{14:.4f}\t{15:.4f}\t{16:.4f}\t{17:.4f}\t{18:.4f}\t{19:.4f}\t{20:.4f}\t{21:.4f}\t{22:.4f}\t{23:.4f}\t{24:.4f}\t{25:.4f}\t{26:.4f}\t{27:.4f}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{28}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{29}\t{30}\t{31}\t{32}\t{33}\t{34}\t{35}\n'.format(i,
                        overstory,rarc,rmin,lai[0],lai[1],lai[2],lai[3],lai[4],lai[5],lai[6],lai[7],lai[8],lai[9],lai[10],lai[11],alb[0],alb[1],alb[2],alb[3],alb[4],alb[5],alb[6],alb[7],alb[8],alb[9],alb[10],alb[11],rough,dis,wind_h,rgl,rad_atten,wind_atten,trunk_ratio,comment))

    except IOError: # except raise an error when it doesn't work
        raise IOError('Cannot write output file, error with output veg library file path') 
        
    return

def main():
    
    n_args = len(sys.argv)
    
    # Check user inputs
    if n_args != 6:
        print "Wrong user input"
        print "Script writes the vegetation library file for the VIC model"
        print "usage: python make_veg_lib.py <land cover raster> <LAI data folder> <albedo data folder> <output veg lib file> <land cover classification scheme>"
        print "Exiting system..."
        sys.exit()
        
    else:
        # check that inputs are correct
        if sys.argv[2][-1] != '/':
            print "LAI DIR SHOULD CONTAIN TRAILING '/'"
            print "fixing it for you..."
            sys.argv[2] = sys.argv[2] + "/"
            
        if sys.argv[3][-1] != '/':
            print "ALBEDO DIR SHOULD CONTAIN TRAILING '/'"
            print "fixing it for you..."
            sys.argv[3] = sys.argv[3] + "/"
        
        # do the process
        make_veg_lib(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    
    return

# Execute the main level program if run as standalone    
if __name__ == "__main__":
    main()
