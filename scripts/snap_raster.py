import sys
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

def snap_raster(inputRas,outputRas,templateRas,subGrid,resample):
    
    # Source
    src = gdal.Open(inputRas,GA_ReadOnly)
    src_proj = src.GetProjection()
    srcXSize = src.RasterXSize
    srcYSize = src.RasterYSize
    
    NoData = int(src.GetRasterBand(1).GetNoDataValue())
    
    dtype = gdal.GetDataTypeName(src.GetRasterBand(1).DataType)
    
    datatypes = {'Byte':gdal.GDT_Byte,'Int16':gdal.GDT_Int16,
                 'UInt16':gdal.GDT_UInt16,'UInt32':gdal.GDT_UInt32,
                 'Int32':gdal.GDT_Int32,'Float32':gdal.GDT_Float32,
                 'Float64':gdal.GDT_Float64}
    
    src_dtype = datatypes[dtype]
    
    # We want a section of source that matches this:
    match_ds = gdal.Open(templateRas,GA_ReadOnly)
    match_proj = match_ds.GetProjection()
    match_geotrans = match_ds.GetGeoTransform()
    matchXSize = match_ds.RasterXSize
    matchYSize = match_ds.RasterYSize
    
    if subGrid == "True":
        xRatio = int(np.round(srcXSize / matchXSize))
        yRatio = int(np.round(srcYSize / matchYSize))
        
        wide = matchXSize * xRatio
        high = matchYSize * yRatio
        
        outGeom = [match_geotrans[0],(match_geotrans[1]/float(xRatio)) ,0,
                   match_geotrans[3],0,(match_geotrans[5]/float(yRatio))]
        
    else:
        wide = matchXSize
        high = matchYSize
        outGeom = match_geotrans
        
    if resample != None:
        samps = {'nearest':GRA_NearestNeighbour,'bilinear':GRA_Bilinear,'cubic':GRA_Cubic,
                 'spline':GRA_CubicSpline, 'mean':GRA_Average, 'mode':GRA_Mode}	
        sampMethod = samps[resample]
    
    # Output / destination
    dst = gdal.GetDriverByName('GTiff').Create(outputRas, wide, high, 1, src_dtype)
    dst.SetGeoTransform(outGeom)
    dst.SetProjection(match_proj)
    band = dst.GetRasterBand(1)
    band.SetNoDataValue(NoData)
    
    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, sampMethod)
    
    dst = None # Flush
    band = None
    
    return

def main():
    
    snap_raster(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])

    return


if __name__ == "__main__":
    main()
