#******************************************************************************
# FILE: calibrate_vic.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 22 Feb. 2017
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script performs a simple calibration process by inserting random
#          parameterization into the VIC model  
# DEPENDENCIES: numpy, pandas, scipy, osgeo (gdal)
#******************************************************************************

# import dependencies
import os
import sys
import datetime
import numpy as np
import pandas as pd
from scipy import stats
from format_soil_params import *
from rout_vic import *
from flux2nc import *


def calibrate_vic(niter):
    
    # define script file path for relative path definitions
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    os.chdir(__location__)
    
    # output calibration table
    outCal = '../data/output/calibration_table.csv'
    
    # specify file paths to pass into system commands
    globalFile = os.path.join(__location__,'../data/input/global.params')
    soilFile = os.path.join(__location__,'../data/input/soil.param')
    gridRas = os.path.join(__location__,'../data/input/gis/Nyando_grid.tif')
    elvRas = os.path.join(__location__,'../data/input/gis/Nyando_basin_ElvAvg.tif')
    slopeRas = os.path.join(__location__,'../data/input/gis/Nyando_basin_SlopeAvg.tif')
    precipRas = os.path.join(__location__,'../data/input/gis/Nyando_basin_PrecipSnap.tif')
    soilRas =os.path.join(__location__,'../data/input/gis/Nyando_basin_SoilAgg.tif')
    
    #gridding paths
    fluxPath = '../data/output/fluxes/'
    fluxNc = '../data/output/'
    
    # specify variables for routing
    fracRas = os.path.join(__location__,'../data/input/gis/Nyando_fraction.tif')
    uhFile = '../data/input//Nyando_UHFile.csv'
    roFile = '../data/output/runoff_2005.nc'
    bfFile = '../data/output/base_2005.nc'
    routOut = '../data/output/rout_out.csv'
    calStart = '20050101'
    calEnd = '20091231'
    
    # specify validation time series path
    valFile = '../data/input/Nyando_discharge.xlsx'    

    # specify range of acceptable parameter
    b_range = [0, 0.5]
    Ws_range = [0,1]
    Ds_range = [0,1]
    s_range = [0.3,1.5]
    
    #empty lists to add calibration results
    b_vals = []
    Ws_vals = []
    Ds_vals = []
    s2_vals = []
    s3_vals = []
    NSEs = []
    Rs = []
    Biases = []
    RMSEs = []
    
    # calibration start time
    t1 = datetime.datetime.now()
    
    # run n number of iterationa
    for i in range(int(niter)):
        
        print "Iteration {0} of {1}".format(i+1,niter) # 
        
        # get random variable for each parameter
        b_val = np.random.uniform(b_range[0],b_range[1],1)[0]
        Ws_val = np.random.uniform(Ws_range[0],Ws_range[1],1)[0]
        Ds_val = np.random.uniform(Ds_range[0],Ds_range[1],1)[0]
        s2_val = np.random.uniform(s_range[0],s_range[1],1)[0]
        s3_val = np.random.uniform(s_range[0],s_range[1],1)[0]
 
        # use random parameters in soil parameter file
        format_soil_params(gridRas,soilRas,elvRas,precipRas,slopeRas,soilFile,
                           b_val, Ws_val, Ds_val, s2_val, s3_val)
                           
        # run the VIC model
        os.system('../data/vicNl -g {0} 2> ../data/output/vic.log'.format(globalFile))
        
        # grid flux
        flux2nc(fluxPath,fluxNc, 5,2005,2009)
        flux2nc(fluxPath,fluxNc, 6,2005,2009)
        
        # rout the VIC model
        rout_vic(uhFile,fracRas,roFile,bfFile,routOut,calStart,calEnd,daily='False')
        
        # read simulated and observed data for time period
        simCsv = pd.read_csv(routOut)
        obsData = pd.ExcelFile(valFile)
        
        obsSheet = obsData.parse('Monthly')
        
        obsSeries = np.array(obsSheet.Q)[4:59]
        
        simSeries = np.array(simCsv.Discharge)[5:]
        
        # calculate model performance statistics
        r = stats.pearsonr(obsSeries,simSeries)
        nse = 1 - (sum((obsSeries-simSeries)**2)/sum((obsSeries-obsSeries.mean())**2))
        bias = np.mean(simSeries-obsSeries)
        rmse = np.mean(np.abs(obsSeries-simSeries))
        
        # append variables to lists
        b_vals.append(b_val)
        Ws_vals.append(Ws_val)
        Ds_vals.append(Ds_val)
        s2_vals.append(s2_val)
        s3_vals.append(s3_val)
        NSEs.append(nse)
        Rs.append(r[0])
        Biases.append(bias)
        RMSEs.append(rmse)
        
    deltat = datetime.datetime.now()-t1
    print 'Processing time for {0} iterations: {1}'.format(niter,deltat)
        
    # create dictionary to be put in a dataframe
    d = {'Infilt':b_vals,"Ws":Ws_vals,"Ds":Ds_vals,"S2Depth":s2_vals,"S3Depth":s3_vals,
         "NSE":NSEs,'R':Rs,"Bias":Biases,'RMSE':RMSEs}
            
    # convert dictionary to datafram
    df = pd.DataFrame(data=d)
        
    # save dataframe to csv
    df.to_csv(outCal)
    
    return
        
def main():
    
    n_args = len(sys.argv)
    
    # Check user inputs
    if n_args != 1:
        print "Wrong user input"
        print "Script used to perform a simple random calibration process for the VIC model"
        print "usage: python calibrate_vic.py <number of iterations> "
        print "Exiting system..."
        sys.exit()
        
    else: # do the process
        calibrate_vic(sys.argv[1])
        
    return
    
# Execute the main level program if run as standalone    
if __name__ == "__main__":
    main()
        