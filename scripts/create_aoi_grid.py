import os
import sys
import numpy as np
from osgeo import gdal,ogr,osr


def create_aoi_grid(inputAOI,outputGrid,gridSize):
    
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
    
    hiResRatio = 50.
    
    # Create high res source for boundary accuracy
    highResGrid = gridSize / hiResRatio
    
    x_hres = int(np.ceil((x_max - x_min) / highResGrid))
    y_hres = int(np.ceil((y_max - y_min) / highResGrid))
    
    mem_ds = gdal.GetDriverByName('MEM').Create('', x_hres, y_hres, gdal.GDT_Byte)
    mem_ds.SetGeoTransform((x_min, highResGrid, 0, y_max, 0, -highResGrid))
    band = mem_ds.GetRasterBand(1)
    band.SetNoDataValue(NoData_value)
    
    # Rasterize shapefile to high resolution grid
    gdal.RasterizeLayer(mem_ds, [1], source_layer, burn_values=[1])
    
    array = band.ReadAsArray()
    
    mem_ds = None
    band = None
    
    # Create the destination data source
    x_res = int(np.ceil((x_max - x_min) / gridSize))
    y_res = int(np.ceil((y_max - y_min) / gridSize))
    drv =  gdal.GetDriverByName('GTiff')
    target_ds = drv.Create(os.path.join(__location__,outputGrid), x_res, y_res, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((x_min, gridSize, 0, y_max, 0, -gridSize))
    target_ds.SetProjection(spatialRef.ExportToWkt())
    
    outMask = np.zeros([y_res,x_res])
    
    for i in range(y_res):
        y1 = int(i*hiResRatio)
        y2 = int(y1+hiResRatio)
        for j in range(x_res):
            x1 = int(j*hiResRatio)
            x2 = int(x1+hiResRatio)
            
            tmp = array[y1:y2,x1:x2]
            
            if np.any(tmp==1)==True:
                outMask[i,j] = 1
    
    band = target_ds.GetRasterBand(1)
    band.WriteArray(outMask)
    band.SetNoDataValue(NoData_value)
    
    target_ds = None
    band = None
    
    return

def main():
    create_aoi_grid(sys.argv[1],sys.argv[2],float(sys.argv[3]))

    return

if __name__ == "__main__":
    main()
