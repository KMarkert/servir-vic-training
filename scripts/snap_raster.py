#******************************************************************************
# FILE: snap_raster.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 02 Apr. 2017
# PURPOSE: This script aligns the grid cells of two rasters and does resampling
#          if needed
# DEPENDENCIES: numpy, osgeo (gdal)
#******************************************************************************

# import dependencies
import sys
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

def snap_raster(inputRas,outputRas,templateRas,subGrid,resample):
    """
    FUNCTION: snap_raster
    ARGUMENTS: inputRas - input raster file to be snapped with TIF file extension
               outputRas - output snapped raster file with TIF file extension
               templateRas - input raster that the raster will be snapped to with TIF file extension
               subGrid - boolean value to set whether the output raster's 
                         resolution will be at the template raster resoltion or not
               resample - resampling method to be used when snapping
    KEYWORDS: n/a
    RETURNS: n/a
    NOTES: Returns no variables but writes an output file
    """
    
    # try to do the process
    try:
        
        # read source raster
        src = gdal.Open(inputRas,GA_ReadOnly)
        src_proj = src.GetProjection()
        srcXSize = src.RasterXSize
        srcYSize = src.RasterYSize
        
        # get no data value
        NoData = int(src.GetRasterBand(1).GetNoDataValue())
        
        # get source raster datatype
        dtype = gdal.GetDataTypeName(src.GetRasterBand(1).DataType)
        
        # get destination data type through a lookup dictionary
        datatypes = {'Byte':gdal.GDT_Byte,'Int16':gdal.GDT_Int16,
                    'UInt16':gdal.GDT_UInt16,'UInt32':gdal.GDT_UInt32,
                    'Int32':gdal.GDT_Int32,'Float32':gdal.GDT_Float32,
                    'Float64':gdal.GDT_Float64}
                    
        src_dtype = datatypes[dtype]
        
        # try getting the resampling method        
        try:
            samps = {'nearest':GRA_NearestNeighbour,'bilinear':GRA_Bilinear,'cubic':GRA_Cubic,
                    'spline':GRA_CubicSpline, 'mean':GRA_Average, 'mode':GRA_Mode}	
            sampMethod = samps[resample]
        # if the resampling method does not exist, give an error
        except KeyError:
            raise KeyError('{0} is not a valid resampling method'.format(resample))
            
        # check to make sure that the output data type will make sense with the resampling method
        if (resample != 'nearest') and (resample != 'mode'):
            if 'Int' in dtype:
                src_dtype = datatypes['Float32']
            elif 'Byte' in dtype:
                src_dtype = datatypes['Float32']
            else:
                pass
    
        # We want a section of source that matches this:
        match_ds = gdal.Open(templateRas,GA_ReadOnly)
        match_proj = match_ds.GetProjection()
        match_geotrans = match_ds.GetGeoTransform()
        matchXSize = match_ds.RasterXSize
        matchYSize = match_ds.RasterYSize
        
        # check if the subGrid argument is true or no
        if (subGrid == "True") | (subGrid == "true") | (subGrid == "TRUE"):
            # calculate the new raster resolution closest to source resolution
            xRatio = int(np.round(srcXSize / matchXSize))
            yRatio = int(np.round(srcYSize / matchYSize))
            
            wide = matchXSize * xRatio
            high = matchYSize * yRatio
            
            # create ouput geotransform with new raster resolution
            outGeom = [match_geotrans[0],(match_geotrans[1]/float(xRatio)) ,0,
                    match_geotrans[3],0,(match_geotrans[5]/float(yRatio))]
            
        else:
            # set output geotransform to the template raster
            wide = matchXSize
            high = matchYSize
            outGeom = match_geotrans
        
        # create output / destination raster
        dst = gdal.GetDriverByName('GTiff').Create(outputRas, wide, high, 1, src_dtype)
        dst.SetGeoTransform(outGeom)
        dst.SetProjection(match_proj)
        band = dst.GetRasterBand(1)
        band.SetNoDataValue(NoData)
        
        # Do the work
        gdal.ReprojectImage(src, dst, src_proj, match_proj, sampMethod)
        
        # flush the dataset
        dst = None 
        band = None
    
    # if not working, give error message
    except AttributeError:
        raise IOError('Raster file input error, check that all paths are correct')
    
    return

# main level program
def main():
    
    n_args = len(sys.argv)
    
    # Check user inputs
    if n_args != 6:
        print "Wrong user input"
        print "Script converts shapefile to raster grid"
        print "usage: python snap_raster.py <input raster> <output raster> <template raster> <sub grid boolean> <resampling method>"
        print "Exiting system..."
        sys.exit()
    
    else:
        # Pass command line arguments into function
        snap_raster(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])

    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()
