#******************************************************************************
# FILE: format_snow_params.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 03 Apr. 2017
# PURPOSEe This script takes geotiff data for snow band parameterization and converts
#          it to the needed text format for the VIC model
# DEPENDENCIES: numpy, osgeo (gdal)
#******************************************************************************

# import dependencies
from __future__ import print_function
import os
import sys
import warnings
import numpy as np
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *

# set system to ignore simple warnings
warnings.simplefilter("ignore")

def format_snow_params(basinMask,elvHiRes,outSnow,interval):
    """
    FUNCTION: format_snow_params
    ARGUMENTS: basinMask - path to template raster to run VIC model at
               elvHiRes - path elevation raster dataset at native resolution
               outsnow - path output snow parameter file
               interval - vertical distance to do equal interval segmentation
    KEYWORDS: n/a
    RETURNS: n/a
    NOTES: Does not return a variable but writes an output file
    """

    band = 1 # constant variable for reading in data

    interval = int(interval) # force equal interval value to be int type

    # make a list of input raster files
    infiles = [basinMask,elvHiRes]

    # try to read in the raster data
    try:

        # read basin grid raster
        ds = gdal.Open(infiles[0],GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        mask = BandReadAsArray(b1)
        maskRes = ds.GetGeoTransform()[1]
        ds = None
        b1 = None

        # read hi res elevation raster
        ds = gdal.Open(infiles[1],GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        elvhires = BandReadAsArray(b1)
        clsRes = ds.GetGeoTransform()[1]
        ds = None
        b1 = None

    # if not working, give error message
    except AttributeError:
        raise IOError('Raster file input error, check that all paths are correct')

    # mask elevation values less than 0
    elvhires[np.where(elvhires<0)] = np.nan

    # get ratio of high resoltion to low resolution
    clsRatio = int(maskRes/clsRes)

    # check if the output parameter file exists, if so delete it
    if os.path.exists(outSnow)==True:
        os.remove(outSnow)

    nbands = [] # blank list

    # try to write to output snow parameter file
    try:

        with open(outSnow, 'w') as f:

            cnt = 1 # set grid cell id counter

            # pass counter
            pass_counter = range(2)

            # perform two passes on the raster data
            # 1) to grab the maximum number of bands for a given pixel
            # 2) to calculate the snow band parameters and write to output file
            for pass_cnt in pass_counter:

                # loop over each pixel in the template raster
                for i in range(mask.shape[0]):
                    cy1 = i*clsRatio
                    cy2 = cy1+clsRatio

                    for j in range(mask.shape[1]):
                        cx1 = j*clsRatio
                        cx2 = cx1+clsRatio

                        # get all hi res pixels in a template pixel
                        tmp = elvhires[cy1:cy2,cx1:cx2]
                        if tmp.size == 0:
			    tmp=tmp2[:]
      		        
                        
                        # create blank array for number of bands calculation...
                        if mask[i,j] == 1: # ...if it is not a masked pixel
			    tmp2=tmp

			    if np.all(tmp == np.nan) == True:
				tmp[:,:] = 0
				
                            # find min and max values for interval
                            minelv = np.nanmin(tmp.astype(int)) - (np.nanmin(tmp.astype(int))%interval)
      		            maxelv = np.nanmax(tmp.astype(int)) + (np.nanmax(tmp.astype(int))%interval)
                        
                            #print(np.min(tmp).mask)
                            # create an array of band limits
      		            bands = np.arange(minelv, maxelv+interval,interval)

      		            bcls = np.zeros_like(tmp)
      		            bcls[:,:] = -1

                            # get the number of bands per pixel
          		    for b in range(bands.size-1):

                                bcls[np.where((tmp>=bands[b])&(tmp<bands[b+1]))] = b # band counter

                                # if it's the first pass get number of bands for each pixel
                                if pass_cnt == 0:
                                    uniqcnt = np.unique(bcls[np.where(tmp>0)])
          		            nbands.append(uniqcnt.size) # save to a list for second pass

          		    if pass_cnt == 1:
                                uniqcnt = np.unique(bcls[np.where(tmp>0)])
                                #clscnt = np.bincount(tmp.ravel())

                                f.write('{0}\t'.format(cnt)) # write grid cell id

                                # find frational area for each band and write to file
                                for c in range(maxbands):
                                    try:
                                        idx = np.where(bcls==uniqcnt[c])
                                        num = np.float(bcls[np.where(bcls>=0)].size)
					if num == 0:
					    num = np.float(idx[0].size)
                                        frac = np.float(idx[0].size) / num
                                    except IndexError:
                                        frac = 0

             			    f.write('{0:.4f}\t'.format(frac))

                                # calculate the mean elevation for each band and write to file
                                for c in range(maxbands):
                                    try:
                                        idx = np.where(bcls==uniqcnt[c])
                                        muelv = np.nanmean(tmp[idx])

                                    except IndexError:
                                        muelv = 0

                                    f.write('{0:.4f}\t'.format(muelv))

                                # calculate the precipitation fractions and write to file
      		                for c in range(maxbands):
     			            try:
                                        idx = np.where(bcls==uniqcnt[c])
					num = np.float(bcls[np.where(bcls>=0)].size)
					if num == 0:
					    num = np.float(idx[0].size)
                                        frac = np.float(idx[0].size) / num
                                    except IndexError:
                                        frac = 0

     			            f.write('{0:.4f}\t'.format(frac))

     			        f.write('\n') # write return value for new line

       	                if pass_cnt == 1 & mask[i,j] == 1:
                            cnt += 1 # plus one to the grid cell id counter

                if pass_cnt == 0:
       	            maxbands = max(nbands) # maximum number of bands for a pixel

        # print the number of bands for user to input into global parameter file
        print('Number of maximum bands: {0}'.format(maxbands))

    # except raise an error when it doesn't work
    except IOError:
        raise IOError('Cannot write output file, error with output snow parameter file path')

    return

def main():
    n_args = len(sys.argv)

    # Check user inputs
    if n_args != 5:
        print ("Wrong user input")
        print ("Script writes the snow band parameter file for the VIC model")
        print ("usage: python format_snow_params.py <template raster> <elevation raster> <output snow band file> <interval for snow bands>")
        print ("Exiting system...")
        sys.exit()

    else:
        # Pass command line arguments into function
        format_snow_params(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()
