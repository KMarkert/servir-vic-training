#******************************************************************************
# FILE: format_soil_params.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 02 Apr. 2017
# PURPOSE: This script takes geotiff data for soil parameterization and converts
#          it to the needed text format for the VIC model
# DEPENDENCIES: numpy, pandas, osgeo (gdal)
#******************************************************************************

from __future__ import print_function
import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *

# set system to ignore simple warnings
warnings.simplefilter("ignore")

def get_soil_params(scls, sdata,subsoil):
    """
    FUNCTION: get_soil_params
    ARGUMENTS: scls - soil class value from HWSD raster file
               sdata - HWSD data table for the top layer
               subsoil - HWSD data table for the bottom layer
    KEYWORDS: n/a
    RETURNS: dictionary with soil attribute information
    NOTES: n/a
    """

    # grab an array of potential usda soil classes
    usdaclss = sdata[1][np.where(sdata[0]==scls)]
    susdaclss = subsoil[1][np.where(subsoil[0]==scls)]

    # find the modal usda class information
    try: usdacls = np.bincount(usdaclss).argmax()
    # if there is an error, set a default value
    except ValueError: usdacls = 9

    # sub soil usda class
    try: susdacls = np.bincount(susdaclss).argmax()
    except ValueError: susdacls = 9

    # get bulk density information
    tbden = np.nanmean(sdata[2][np.where(sdata[0]==scls)]) * 1000.
    sbden = np.nanmean(subsoil[2][np.where(subsoil[0]==scls)]) * 1000.

    # check if bulk density information is true, if not, then set default value
    if np.isnan(tbden)==True:
        tbden = 1331.
    if np.isnan(sbden)==True:
        sbden = 1395.

    # get organic content information
    toc = np.nanmean(sdata[3][np.where(sdata[0]==scls)]) / 100.
    soc = np.nanmean(subsoil[3][np.where(subsoil[0]==scls)]) / 100.

    # check if organic content information is true, if not, then set default value
    if np.isnan(toc)==True:
        tbden = 0.020
    if np.isnan(soc)==True:
        soc = 0.015

    # grab an array of drainage values
    drn = sdata[4][np.where(sdata[0]==scls)]

    # reduce to modal drainage value
    try: drncls = np.bincount(drn).argmax()
    except ValueError: drncls = 4

    # return dictionary of soil attributes
    return {'topUSDA':usdacls,'subUSDA':susdacls,'topBulkDen':tbden,
            'subBulkDen':sbden,'topOC':toc,'subOC':soc,'drainage':drncls}

def format_soil_params(basinMask,HWSD,basinElv,AnnPrecip,Slope,outsoil,
                       b_val=None,Ws_val=None,Ds_val=None,s2=None,s3=None):

    band = 1 # constant variable for reading in data

    # define script file path for relative path definitions
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # define soil type lookup file path
    attriFile = os.path.join(__location__,'soil_type_attributes.json')

    # open/read soil type json file
    with open(attriFile) as data_file:
        attriData = json.load(data_file)

    # pass look up information into variable
    soilAttributes = attriData['classAttributes']

    # define drainage type lookup file path
    attriFile = os.path.join(__location__,'drain_type_attributes.json')

    # open/read drainage type json file
    with open(attriFile) as data_file:
        attriData = json.load(data_file)

    # pass lookup information into variable
    drainAttributes = attriData['classAttributes']

    # define path to HWSD table
    csvfile = os.path.join(__location__,'HWSD_CLS_DATA.csv')

    # open and read data for...
    indata = pd.read_csv(csvfile)

    # ...top soil layer...
    soildata = [np.array(indata.MU_GLOBAL),
                np.array(indata.T_USDA_TEX_CLASS,dtype=np.int32),
                np.array(indata.T_BULK_DENSITY,dtype=np.float),
                np.array(indata.T_OC,dtype=np.float),
                np.array(indata.DRAINAGE,dtype=np.int32)]
    # ...and bottom soil layer
    subsoil = [np.array(indata.MU_GLOBAL),
               np.array(indata.S_USDA_TEX_CLASS,dtype=np.int32),
               np.array(indata.S_BULK_DENSITY,dtype=np.float),
               np.array(indata.S_OC,dtype=np.float)]

    # create list of input raster files
    infiles = [os.path.join(__location__,basinMask),
               os.path.join(__location__,HWSD),
               os.path.join(__location__,basinElv),
               os.path.join(__location__,AnnPrecip),
               os.path.join(__location__,Slope)]

    try:
        # loop over each raster file and pass into an array
        for i in range(len(infiles)):

            # open and read file
            ds = gdal.Open(infiles[i],GA_ReadOnly)
            b1 = ds.GetRasterBand(band)
            var = BandReadAsArray(b1)

            # if it is the first iteration, define blank data array and grab geographic information
            if i == 0:
                data = np.zeros((ds.RasterYSize,ds.RasterXSize, len(infiles)))
                # get geotransform
                gt = ds.GetGeoTransform()
                # define raster bounding box
                lon0 = gt[0] + (gt[1] / 2.)
                lon1 = gt[0] + (data.shape[1]*gt[1])
                lat0 = gt[3] + (data.shape[0]*gt[-1])
                lat1 = gt[3] + (gt[-1] /2)

            # get the no data value
            if i == 2:
                NoData = b1.GetNoDataValue()

            # add raster data to array
            data[:,:,i] = var[:,:]

            # flush variables
            ds = None
            b1 = None

    # if not working, give error message
    except AttributeError:
        raise IOError('Raster file input error, check that all paths are correct')

    # create arrays for all lat/lon points
    lons = np.linspace(lon0,lon1,data.shape[1])
    lats = np.linspace(lat0,lat1,data.shape[0])

    # mesh the lat/lon points
    xx,yy = np.meshgrid(lons,lats)
    yy = np.flipud(yy) # invert the lat

    # define the path to the output soil parameter file
    soilfile = os.path.join(__location__,outsoil)

    # if soil parameter file exists, then delete the file
    if os.path.exists(soilfile)==True:
        os.remove(soilfile)

    cells = 0 # counter

    # try to write the output parameter file
    try:

        # open soil parameter file for writing
        with open(soilfile, 'w') as f:

            cnt = 1 # grid cell id variable

            # loop over each lat/lon point
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):

                    # check if there is data at the lat/lon
                    run = int(data[i,j,0])
                    if run <= 0:
                        run=0
                    if data[i,j,2] == NoData:
                        run = 0

                    # check if to write data for the file based on data
                    if run == 0:
                        # if there is no data, add 1 to cell id variable
                        cnt+=1

                    else:
                        # get soil class and attributes
                        hwsdcls = data[i,j,1]
                        soildic = get_soil_params(hwsdcls, soildata,subsoil)
                        cells+=1

                        # extract data tree from lookup json files
                        soilDrain = drainAttributes[soildic['drainage']-1]['properties']
                        topSoilPro = soilAttributes[soildic['topUSDA']-1]['properties']
                        subSoilPro = soilAttributes[soildic['subUSDA']-1]['properties']

                        # if keywords are not set then pass data from lookup json files
                        if b_val == None:
                            b_val = soilDrain['infilt']
                        if Ds_val == None:
                            Ds_val = soilDrain['Ds']
                        if Ws_val == None:
                            Ws_val = soilDrain['Ws']
                        if s2 == None:
                            s2 = 1.50
                        if s3 == None:
                            s3 = 0.30

                        # start passing information to simple variables
                        grdc = cnt # grid cell id
                        lat = yy[i,j] # latitude
                        lon = xx[i,j] # longitude
                        infilt = b_val # variable infiltration curve parameter
                        Ds = Ds_val # Ds value
                        Dsmax = (data[i,j,4]/100.) * (float(subSoilPro['SatHydraulicCapacity'])*240) # Dsmax value
                        Ws = Ws_val # Ws value
                        c = 2 # exponent in baseflow curve
                        expt = 3+(2*float(topSoilPro['SlopeRCurve'])) # top layer exponent value
                        expt1 = 3+(2*float(subSoilPro['SlopeRCurve'])) # bottom layer exponent value
                        tksat = (float(topSoilPro['SatHydraulicCapacity'])*240) # top layer Ksat value
                        sksat = (float(subSoilPro['SatHydraulicCapacity'])*240) # bottom layer Ksat value
                        phis = -999 # fill value
                        elev = data[i,j,2] # average elevation of gridcell
                        depth = 0.10 # top layer soil depth
                        depth1 =s2 # second layer soil depth
                        depth2 =s3 # bottom layer soil depth
                        avg_t = 27 # average temperature of soil
                        dp = 4 # depth that soil temp does not change
                        tbub = topSoilPro['BubblingPressure'] # top layer bubbling pressure
                        sbub = subSoilPro['BubblingPressure'] # bottom layer bubbling pressure
                        quartz = topSoilPro['Quartz'] # top layer percent quartz
                        quartz1 = subSoilPro['Quartz'] # bottom layer percent quartz
                        bulk_den = float(soildic['topBulkDen']) # top layer bulk density
                        bulk_den1 =float(soildic['subBulkDen']) # bottom layer bulk density
                        soil_den = 2650. # top layer soil density
                        soil_den1 = 2685. # bottom layer soil density
                        t_oc = float(soildic['topOC']) # top layer organic content
                        s_oc = float(soildic['subOC']) # bottom layer organic content
                        org_bulk_den = 0.25*bulk_den # top layer organic bul density
                        org_bulk_den1 = 0.25*bulk_den1 # bottom layer organic bulk density
                        org_soil_den = 1295. # top layer organic soil density
                        org_soil_den1 = 1300. # bottom layer organic soil density
                        off_gmt = lon * 24 / 360. # time zone offset from GMT
                        wrc_frac = (float(topSoilPro['FieldCapacity'])/
                                    float(topSoilPro['Porosity'])) # top layer critical point
                        wrc_frac1 = (float(subSoilPro['FieldCapacity'])/
                                    float(subSoilPro['Porosity'])) # bottom layer critical point
                        wpwp_frac = (float(topSoilPro['WiltingPoint'])/
                                    float(topSoilPro['Porosity'])) # top layer wilting point
                        wpwp_frac1 = (float(subSoilPro['WiltingPoint'])/
                                    float(subSoilPro['Porosity'])) # bottom layer wilting point
                        rough = 0.01 # bare soil roughness coefficient
                        srough = 0.001 # snow roughness coefficient
                        annprecip = data[i,j,3] # climotological average precipitation
                        resid = topSoilPro['Residual'] # top layer residual moisture
                        resid1 = subSoilPro['Residual'] # bottom layer residual moisture
                        fs_act = 1 # boolean value to run frozen soil algorithm
                        init_moist = (bulk_den / soil_den) * depth *1000 # top layer inital moisture conditions
                        initmoist2 = (bulk_den1 / soil_den1) * depth1 *1000 # second layer initial moisture conditions
                        initmoist3 = (bulk_den1 / soil_den1) * depth2 *1000 # bottom layer initial moisture conditions

                        # write the soil parameterization information for each grid cell as a line
                        f.write('{0}\t{1}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8}\t{9}\t{10}\t{10}\t{11}\t{12}\t{12}\t{13}\t{13}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{24}\t{25}\t{26}\t{26}\t{27}\t{28}\t{28}\t{29}\t{30}\t{30}\t{31}\t{32}\t{32}\t{33}\t{34}\t{34}\t{35}\t{36}\t{36}\t{37}\t{38}\t{39}\t{39}\t{40}\t{41}\t{41}\t{42}\t{43}\t{44}\t{45}\t{46}\t{46}\t{47}\n'.format(run,
                                            grdc,lat,lon,infilt,Ds,Dsmax,Ws,c,expt,expt1,tksat,sksat,
                                            phis,init_moist,initmoist2,initmoist3,elev,depth,
                                            depth1,depth2,avg_t,dp,tbub,sbub,quartz,quartz1,bulk_den,
                                            bulk_den1,soil_den,soil_den1,t_oc,s_oc,org_bulk_den,org_bulk_den1,
                                            org_soil_den,org_soil_den1,off_gmt,wrc_frac,wrc_frac1,wpwp_frac,
                                            wpwp_frac1,rough,srough,annprecip,resid,resid1,fs_act))

                        cnt+=1 # plus one to the grid cell id

    # except raise an error
    except IOError:
        raise IOError('Cannot write output file, error with output soil parameter file path')

    return

def main():

    n_args = len(sys.argv)

    # Check user inputs
    if n_args != 7:
        print("Wrong user input")
        print("Script writes the soil parameter file for the VIC model")
        print("usage: python format_soil_params.py <template raster> <soil classification raster> <elevation raster> <annual precip raster> <slope raster> <ouput soil parameter file>")
        print("Exiting system...")
        sys.exit()

    else:
        # Pass command line arguments into function
        format_soil_params(sys.argv[1],sys.argv[2],sys.argv[3],
                           sys.argv[4],sys.argv[5],sys.argv[6])

    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()
