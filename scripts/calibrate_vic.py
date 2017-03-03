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
    
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    os.chdir(__location__)
    
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
    
    for i in range(int(niter)):
        
        print "Iteration {0} of {1}".format(i+1,niter)
        
        b_val = np.random.uniform(b_range[0],b_range[1],1)[0]
        Ws_val = np.random.uniform(Ws_range[0],Ws_range[1],1)[0]
        Ds_val = np.random.uniform(Ds_range[0],Ds_range[1],1)[0]
        s2_val = np.random.uniform(s_range[0],s_range[1],1)[0]
        s3_val = np.random.uniform(s_range[0],s_range[1],1)[0]
        
        #b_val  = 0.06639665    #0.2325047434
        #Ws_val = 0.19130492    #0.0032642209
        #Ds_val = 0.1624679     #0.4410260563
        #s2_val = 0.53565465    #1.4068091083
        #s3_val = 0.30801951    #0.6685944592    
        
        format_soil_params(gridRas,soilRas,elvRas,precipRas,slopeRas,soilFile,
                           b_val, Ws_val, Ds_val, s2_val, s3_val)
        
        os.system('../data/vicNl -g {0} 2> ../data/output/vic.log'.format(globalFile))
        
        flux2nc(fluxPath,fluxNc, 5,2005,2009)
        flux2nc(fluxPath,fluxNc, 6,2005,2009)
        
        rout_vic(uhFile,fracRas,roFile,bfFile,routOut,calStart,calEnd,daily='False')
        
        simCsv = pd.read_csv(routOut)
        
        obsData = pd.ExcelFile(valFile)
        
        obsSheet = obsData.parse('Monthly')
        
        obsSeries = np.array(obsSheet.Q)[4:59]
        #obsSeries = np.array(obsSheet.Q)[120:1795]
        
        simSeries = np.array(simCsv.Discharge)[5:]
        #simSeries = np.array(simCsv.Discharge)[151:]
        
        r = stats.pearsonr(obsSeries,simSeries)
        nse = 1 - (sum((obsSeries-simSeries)**2)/sum((obsSeries-obsSeries.mean())**2))
        bias = np.mean(simSeries-obsSeries)
        rmse = np.mean(np.abs(obsSeries-simSeries))
        
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
        
        
    d = {'Infilt':b_vals,"Ws":Ws_vals,"Ds":Ds_vals,"S2Depth":s2_vals,"S3Depth":s3_vals,
         "NSE":NSEs,'R':Rs,"Bias":Biases,'RMSE':RMSEs}
            
    df = pd.DataFrame(data=d)
        
    df.to_csv(outCal)
    
    return
        
def main():
    calibrate_vic(sys.argv[1])
    return
    
if __name__ == "__main__":
    main()
        