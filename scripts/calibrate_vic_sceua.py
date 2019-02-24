#******************************************************************************
# FILE: calibrate_vic_sceua.py
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

from __future__ import print_function
import os
import sys
import spotpy
import datetime
import numpy as np
import pandas as pd
import xarray as xr
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
        valFile = '../input/Nyando_discharge.xlsx'

        obsData = pd.ExcelFile(valFile)
        obsSheet = obsData.parse('Daily')
        obstimes = pd.date_range('2005-02-01','2013-12-31',freq='D')
        obsSeries = xr.DataArray(obsSheet.Q,coords=[obstimes],dims=['time']).sel(time=slice('2005-03-01','2009-12-31')).data

        self.observations = obsSeries

        return

    def run_vic(self,binfilt=None,Ws=None,Ds=None,c=None,soil_d2=None,soil_d3=None):

        __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__)))

        os.chdir(__location__)

        # specify file paths to pass into system commands
        globalFile = os.path.join(__location__,'../input/global.params')
        soilFile = os.path.join(__location__,'../input/soil.param')
        gridRas = os.path.join(__location__,'../input/gis/Nyando_grid.tif')
        elvRas = os.path.join(__location__,'../input/gis/Nyando_basin_ElvAvg.tif')
        slopeRas = os.path.join(__location__,'../input/gis/Nyando_basin_SlopeAvg.tif')
        precipRas = os.path.join(__location__,'../input/gis/Nyando_basin_Precip.tif')
        soilRas =os.path.join(__location__,'../input/gis/Nyando_basin_SoilsAgg.tif')

        #gridding paths
        fluxPath = '../output/fluxes/calibrated/'
        fluxNc = '../output/netCDFs/calibrated/'

        # specify variables for routing
        fracRas = os.path.join(__location__,'../input/gis/Nyando_gauge_fraction.tif')
        uhFile = '../input/routing/Nyando_UHFile.csv'
        roFile = '../output/netCDFs/calibrated/runoff_2005.nc'
        bfFile = '../output/netCDFs/calibrated/base_2005.nc'
        routOut = '../output/cal_rout_out.csv'
        calStart = self.st.isoformat().replace('-','')
        calEnd = self.et.isoformat().replace('-','')

        # create soil parameter file with new set of variables
        format_soil_params(gridRas,soilRas,elvRas,precipRas,slopeRas,soilFile,
                        binfilt,Ws,Ds,c,soil_d2,soil_d3)

        # execute the vic model
        os.system('../input/vicNl -g {0} 2> ../output/vic.log'.format(globalFile))

        # output data formating
        flux2nc(fluxPath,fluxNc, 5,self.st.year,self.et.year)
        flux2nc(fluxPath,fluxNc, 6,self.st.year,self.et.year)

        # run routing model
        rout_vic(uhFile,fracRas,roFile,bfFile,routOut,calStart,calEnd,daily='True')

        simCsv = pd.read_csv(routOut)
        simtimes = pd.date_range('2005-01-01','2009-12-31',freq='D')
        simSeries = xr.DataArray(simCsv.Discharge,coords=[simtimes],dims=['time'])
        simSeries = simSeries.sel(time=slice('2005-03-01','2009-12-31')).data

        return simSeries


class spotpy_setup(object):
    def __init__(self):
        datastart     = datetime.date(2005,1,1) # calibration start
        dataend       = datetime.date(2009,12,31) # calibration end
        self.vicmodel = vic_model(datastart,dataend,) # routine to run model
        # model parameters to calibrate
        self.params = [spotpy.parameter.Uniform('binfil',0,0.5),
                       spotpy.parameter.Uniform('Ws',0.5,1.),
                       spotpy.parameter.Uniform('Ds',0.,0.5),
                       spotpy.parameter.Uniform('C',0.,5.),
                       spotpy.parameter.Uniform('Soil D2',0.3,1.5),
                       spotpy.parameter.Uniform('Soil D3',0.3,1.5),
                       ]
        return

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self,vector):
        simulations= self.vicmodel.run_vic(binfilt=vector[0],Ws=vector[1],Ds=vector[2],c=vector[3],soil_d2=vector[4],soil_d3=vector[5])
        return simulations

    def evaluation(self,evaldates=False):
        self.vicmodel.get_obs()
        return self.vicmodel.observations

    def objectivefunction(self,simulation,evaluation):
        objectivefunction= -spotpy.objectivefunctions.rmse(evaluation,simulation)
        return objectivefunction

def findBestSim(dbPath):
    csv = pd.read_csv(dbPath+'.csv')

    results = np.array(csv)

    likes = np.array(csv.like1)

    idx = np.abs(likes).argmin()

    b = results[idx,1]
    w = results[idx,2]
    d = results[idx,3]
    c = results[idx,4]
    s2 = results[idx,5]
    s3 = results[idx,6]

    params = [b,w,d,c,s2,s3]

    return params

def runStats(sim,obs):
    nse = 1 - (np.nansum((sim-obs)**2)/np.nansum((obs-obs.mean())**2))
    bias = np.nanmean(sim-obs)
    rmse = np.nanmean(np.sqrt((sim-obs)**2))

    return nse, bias, rmse

def calibrate():
    # calibration setup object
    cal_setup = spotpy_setup()

    outCal = '../output/SCEUA_VIC_Nyando'

    # initialize calibration algorithm with
    sampler = spotpy.algorithms.sceua(cal_setup,dbname=outCal,dbformat='csv')

    results = [] # empty list to append iteration results

    # run calibration process
    sampler.sample(int(sys.argv[1]),ngs=5)
    results.append(sampler.getdata())

    params = findBestSim(outCal)
    print('Best parameter set: {0}'.format(params))

    lastRun = vic_model(datetime.date(2005,1,1),datetime.date(2009,12,31))

    sim = lastRun.run_vic(params[0],params[1],params[2],
                            params[3],params[4],params[5])

    lastRun.get_obs()
    obs = lastRun.observations

    print(sim,obs)
    stats = runStats(sim,obs)
    print('Calibration error statistics...\n \
            NSE:{0}\tBias:{1}\tRMSE{2}'.format(stats[0],stats[1],stats[2]))

    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    calibrate()
