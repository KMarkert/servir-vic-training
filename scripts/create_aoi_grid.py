#******************************************************************************
# FILE: create_aoi_grid.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 02 Apr. 2017
# PURPOSE: This script creates an raster grid over a input shapefile with a user
#          specified cell size
# DEPENDENCIES: numpy, osgeo (gdal,ogr,osr)
#******************************************************************************

# import dependencies
import os
import sys
import numpy as np
from osgeo import gdal,ogr,osr


def create_aoi_grid(inputAOI,outputGrid,gridSize):
    """
    FUNCTION: create_aoi_grid
    ARGUMENTS: inputAOI - input shapefile with SHP file extension
               outputGrid - output raster grid with TIFF file extension
               gridSize - output raster grid cell size in decimal degrees
    KEYWORDS: n/a
    RETURNS: n/a
    NOTES: Returns no variables but writes an output file
    """
    
    # try to do the process
    try:
    
        # define script file path for relative path definitions
        __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__)))
        
        # Define pixel_size and NoData value of new raster
        NoData_value = 0
        
        #Define output coordinate system
        spatialRef = osr.SpatialReference()
        spatialRef.SetWellKnownGeogCS('WGS_84')
        
        # Open the data source and read in the extent
        source_ds = ogr.Open(os.path.join(__location__,inputAOI))
        source_layer = source_ds.GetLayer()
        x_min, x_max, y_min, y_max = source_layer.GetExtent()
        
        # Create high res source for boundary accuracy
        hiResRatio = 50.
        highResGrid = gridSize / hiResRatio
        
        # Get high res cell size
        x_hres = int(np.ceil((x_max - x_min) / highResGrid))
        y_hres = int(np.ceil((y_max - y_min) / highResGrid))
        
        # Create high res raster in memory
        mem_ds = gdal.GetDriverByName('MEM').Create('', x_hres, y_hres, gdal.GDT_Byte)
        mem_ds.SetGeoTransform((x_min, highResGrid, 0, y_max, 0, -highResGrid))
        band = mem_ds.GetRasterBand(1)
        band.SetNoDataValue(NoData_value)
        
        # Rasterize shapefile to high resolution grid
        gdal.RasterizeLayer(mem_ds, [1], source_layer, burn_values=[1])
        
        # Get rasterized high res shapefile
        array = band.ReadAsArray()
        
        # Flush memory file
        mem_ds = None
        band = None
        
        # Create the destination data source
        x_res = int(np.ceil((x_max - x_min) / gridSize))
        y_res = int(np.ceil((y_max - y_min) / gridSize))
        drv =  gdal.GetDriverByName('GTiff')
        target_ds = drv.Create(os.path.join(__location__,outputGrid), x_res, y_res, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform((x_min, gridSize, 0, y_max, 0, -gridSize))
        target_ds.SetProjection(spatialRef.ExportToWkt())
        
        # Create blank mask array at target res
        outMask = np.zeros([y_res,x_res])
        
        # Loop over array to find the elements the high res raster falls in
        for i in range(y_res):
            y1 = int(i*hiResRatio)
            y2 = int(y1+hiResRatio)
            for j in range(x_res):
                x1 = int(j*hiResRatio)
                x2 = int(x1+hiResRatio)
                
                tmp = array[y1:y2,x1:x2]
                
                if np.any(tmp==1)==True:
                    outMask[i,j] = 1
        
        # set the mask array to the target file
        band = target_ds.GetRasterBand(1)
        band.WriteArray(outMask)
        band.SetNoDataValue(NoData_value)
        
        # Flush the target file
        target_ds = None
        band = None
    
    # if not working, then give an error message
    except AttributeError:
        raise IOError('Raster file input error, check that all paths are correct')
    
    return

# Main level program for use in the command line
def main():
    
    n_args = len(sys.argv)
    
    # Check user inputs
    if n_args != 4:
        print "Wrong user input"
        print "Script converts shapefile to raster grid"
        print "usage: python create_aoi_grid.py <input shapefile> <output raster file> <output cell size>"
        print "Exiting system..."
        sys.exit()
    
    else:
        # Pass command line arguments into function
        create_aoi_grid(sys.argv[1],sys.argv[2],float(sys.argv[3]))

    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()
