#******************************************************************************
# FILE: calibrate_vic.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 22 Feb. 2017
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script performs a calibration process for the VIC model using 
#          the Shuffled Complex Evolution (SCE-UA) algorithm
# DEPENDENCIES: numpy, pandas, scipy, spotpy
#******************************************************************************

import os
import sys
import spotpy
import datetime
import numpy as np
import pandas as pd
from format_soil_params import *
from rout_vic import *
from flux2nc import *


class vic_model(object):
    
    def __init__(self,startTime,endTime):
        self.st = startTime
        self.et = endTime
        
        return
        
    def get_obs(self):
        # specify validation time series path
        valFile = '../data/input/Nyando_discharge.xlsx'
        
        obsData = pd.ExcelFile(valFile)
        
        obsSheet = obsData.parse('Monthly')
        
        obsSeries = np.array(obsSheet.Q)[4:59]
        
        self.observations = obsSeries
        
        return

    def run_vic(self,binfilt=None,Ws=None,Ds=None,soil_d2=None,soil_d3=None):
        
        __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__)))
        
        os.chdir(__location__)
                
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
        calStart = self.st.isoformat().replace('-','')
        calEnd = self.et.isoformat().replace('-','')
        
        # create soil parameter file with new set of variables
        format_soil_params(gridRas,soilRas,elvRas,precipRas,slopeRas,soilFile,
                        binfilt,Ws,Ds,soil_d2,soil_d3)
                        
        # execute the vic model
        os.system('../data/vicNl -g {0} 2> ../data/output/vic.log'.format(globalFile))
        
        # output data formating
        flux2nc(fluxPath,fluxNc, 5,self.st.year,self.et.year)
        flux2nc(fluxPath,fluxNc, 6,self.st.year,self.et.year)
        
        # run routing model
        rout_vic(uhFile,fracRas,roFile,bfFile,routOut,calStart,calEnd,daily='False')
        
        simCsv = pd.read_csv(routOut)
        
        simSeries = np.array(simCsv.Discharge)[5:]
            
        return simSeries
        
    
class spotpy_setup(object):
    def __init__(self):
        datastart     = datetime.date(2005,1,1) # calibration start
        dataend       = datetime.date(2009,12,31) # calibration end
        self.vicmodel = vic_model(datastart,dataend,) # routine to run model
        # model parameters to calibrate
        self.params = [spotpy.parameter.Uniform('binfil',0,0.5),
                       spotpy.parameter.Uniform('Ws',0.,1.),
                       spotpy.parameter.Uniform('Ds',0.,1.),
                       spotpy.parameter.Uniform('Soil D2',0.3,1.5),
                       spotpy.parameter.Uniform('Soil D3',0.3,1.5),
                       ]
        return
                       
    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self,vector):
        simulations= self.vicmodel.run_vic(binfilt=vector[0],Ws=vector[1],Ds=vector[2],soil_d2=vector[3],soil_d3=vector[4])
        return simulations
    
    def evaluation(self,evaldates=False):
        self.vicmodel.get_obs()
        return self.vicmodel.observations

    def objectivefunction(self,simulation,evaluation):
        objectivefunction= spotpy.objectivefunctions.nashsutcliff(evaluation,simulation)
        return objectivefunction

        
def calibrate():
    # calibration setup object
    cal_setup = spotpy_setup()
    
    # initialize calibration algorithm with 
    sampler = spotpy.algorithms.sceua(cal_setup,dbname='SCEUA_VIC',dbformat='csv')
    
    results = [] # empty list to append iteration results
    
    # run calibration process
    sampler.sample(int(sys.argv[1]),ngs=5)
    results.append(sampler.getdata)
    
    print spotpy.analyser.get_best_parameterset(results)
    
    evaluation = spotpy_setup().evaluation()
    
    print evaluation
    
    return
    
# Execute the main level program if run as standalone
if __name__ == "__main__":
    calibrate()