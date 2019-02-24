#******************************************************************************
# FILE: format_veg_params.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script takes geotiff data for vegetation parameterization and
#          converts it to the needed text format for the VIC model
# DEPENDENCIES: numpy, osgeo (gdal)
#******************************************************************************

# import dependencies
from __future__ import print_function
import os
import sys
import json
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *

def format_veg_params(basinMask,lcData,outVeg,scheme='IGBP'):
    """
    FUNCTION: format_veg_params
    ARGUMENTS: basinMask - path to basin template raster
               lcData - path to land cover raster
               outveg - path output vegetation parameter file
    KEYWORDS:  scheme - Abbreviation of land cover classification scheme the
                        input land cover data is formatted in
    RETURNS: n/a
    NOTES: Returns no variables but writes an output file
    """

    # define script file path for relative path definitions
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # get land cover classification scheme and create path to lookup table
    if scheme == 'IGBP':
        attriFile = os.path.join(__location__,'veg_type_attributes_igbp.json')
    elif scheme == 'GLCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_glcc.json')
    elif scheme == 'IPCC':
        attriFile = os.path.join(__location__,'veg_type_attributes_ipcc.json')
    else:
        raise SyntaxError('Land cover classification scheme not supported')

    band = 1 # constant variable for reading in data

    # open/read veg scheme json file
    with open(attriFile) as data_file:
        attriData = json.load(data_file)

    # pass look up information into variable
    clsAttributes = attriData['classAttributes']

    # create list of input raster files
    infiles = [os.path.join(__location__,basinMask),
               os.path.join(__location__,lcData)]

    try: # try to read in the raster data

        # read basin grid raster
        ds = gdal.Open(infiles[0],GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        mask = BandReadAsArray(b1)
        maskRes = ds.GetGeoTransform()[1]
        ds = None
        b1 = None

        # read land cover raster
        ds = gdal.Open(infiles[1],GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        lccls = BandReadAsArray(b1)
        clsRes = ds.GetGeoTransform()[1]
        ds = None
        b1 = None

    # if not working, give error message
    except AttributeError:
        raise IOError('Raster file input error, check that all paths are correct')

    ratio = maskRes/clsRes # get ratio of high resoltion to low resolution

    # get file path to output file
    vegfile = os.path.join(__location__,outVeg)

    # check if the output parameter file exists, if so delete it
    if os.path.exists(vegfile)==True:
        os.remove(vegfile)

    try: # try to write output veg parameter file

        # open output file for writing
        with open(vegfile, 'w') as f:

            cnt = 1 # grid cell id counter

            # loop over each pixel in the template raster
            for i in range(mask.shape[0]):
                y1 = int(i*ratio)
                y2 = int(y1+ratio)
                for j in range(mask.shape[1]):
                    x1 = int(j*ratio)
                    x2 = int(x1+ratio)

                    # get land cover data within template raster pixel
                    tmp = lccls[y1:y2,x1:x2]

                    # if not a mask value...
                    if mask[i,j] == 1:
                        # if there are nodata values...
                        if np.any(tmp>len(clsAttributes)-1)==True:
                            # ...find where they are...
                            negdx = np.where(tmp>len(clsAttributes)-1)
                            try:
                                # ...and calculate hisogram for only data values
                                tmp[negdx] = np.bincount(tmp[np.where(tmp<len(
                                                clsAttributes)-1)].ravel()).argmax()
                            except ValueError:
                                tmp[negdx] = 0
                        uniqcnt = np.unique(tmp) # number of unique values
                        clscnt = np.bincount(tmp.ravel())

                        # check if there is only 1 value in unique value array
                        if type(uniqcnt).__name__ == 'int':
                            Nveg = 1 # if so, only 1 value
                        else:
                            Nveg = len(uniqcnt) # if not, n veg equal to unique values

                        # write the grid cell id and number of classes
                        f.write('{0} {1}\n'.format(cnt,Nveg))

                        # if there is vegetation data...
                        if Nveg != 0:
                            # ...loop over each class...
                            for t in range(uniqcnt.size):
                                # ...and grab the parameter attributes
                                vegcls = int(uniqcnt[t]) # class value
                                Cv = np.float(clscnt[uniqcnt[t]])/np.float(clscnt.sum()) # percent coverage
                                attributes = clsAttributes[vegcls]['properties'] # intermediate variale
                                rdepth1 = str(attributes['rootd1']) # rooting depth layer 1
                                rfrac1 = str(attributes['rootfr1']) # rooting fraction layer 1
                                rdepth2 = str(attributes['rootd2']) # rooting depth layer 2
                                rfrac2 = str(attributes['rootfr2']) # rooting fraction layer 2
                                rdepth3 = str(attributes['rootd3']) # rooting depth layer 3
                                rfrac3 = str(attributes['rootfr3']) # rooting fraction layer 3

                                # write the veg class information
                                f.write('\t{0} {1:.4f} {2} {3} {4} {5} {6} {7}\n'.format(vegcls,Cv,rdepth1,rfrac1,rdepth2,rfrac2,rdepth3,rfrac3))

                    cnt+=1 # plus one to grid cell id counter

    # except raise an error when it doesn't work
    except IOError:
        raise IOError('Cannot write output file, error with output veg parameter file path')

    return

def main():
    n_args = len(sys.argv)

    # Check user inputs
    if n_args != 5:
        print("Wrong user input")
        print("Script writes the vegetation parameter file for the VIC model")
        print("usage: python format_veg_params.py <template raster> <land cover raster> <output veg parameter file> <land cover classification scheme>")
        print("Exiting system...")
        sys.exit()

    else:
        # Pass command line arguments into function
        format_veg_params(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()
