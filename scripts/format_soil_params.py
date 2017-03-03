import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
from osgeo import gdal
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *

warnings.simplefilter("ignore")

def get_usda_cls(scls, sdata,subsoil):

    usdaclss = sdata[1][np.where(sdata[0]==scls)]
    susdaclss = subsoil[1][np.where(subsoil[0]==scls)]
    
    try: usdacls = np.bincount(usdaclss).argmax()
    except ValueError: usdacls = 9   
    
    try: susdacls = np.bincount(susdaclss).argmax()
    except ValueError: susdacls = 9
    
    tbden = np.nanmean(sdata[2][np.where(sdata[0]==scls)]) * 1000.
    sbden = np.nanmean(subsoil[2][np.where(subsoil[0]==scls)]) * 1000.
    
    if np.isnan(tbden)==True:
        tbden = 1331.
    if np.isnan(sbden)==True:
        sbden = 1395.
        
    toc = np.nanmean(sdata[3][np.where(sdata[0]==scls)]) / 100.
    soc = np.nanmean(subsoil[3][np.where(subsoil[0]==scls)]) / 100.
    
    if np.isnan(toc)==True:
        tbden = 0.020
    if np.isnan(soc)==True:
        soc = 0.015    
    
    drn = sdata[4][np.where(sdata[0]==scls)]
    
    try: drncls = np.bincount(drn).argmax()
    except ValueError: drncls = 4
        
    return {'topUSDA':usdacls,'subUSDA':susdacls,'topBulkDen':tbden,
            'subBulkDen':sbden,'topOC':toc,'subOC':soc,'drainage':drncls}

def format_soil_params(basinMask,HWSD,basinElv,AnnPrecip,Slope,outsoil,
                       b_val=None,Ws_val=None,Ds_val=None,s2=None,s3=None):
    
    band = 1
    
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    attriFile = os.path.join(__location__,'soil_type_attributes.json')
    
    with open(attriFile) as data_file:    
        attriData = json.load(data_file)
    
    soilAttributes = attriData['classAttributes']
    
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    attriFile = os.path.join(__location__,'drain_type_attributes.json')
    
    with open(attriFile) as data_file:    
        attriData = json.load(data_file)
    
    drainAttributes = attriData['classAttributes']
    
    __location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    csvfile = os.path.join(__location__,'HWSD_CLS_DATA.csv')
    
    indata = pd.read_csv(csvfile)
    
    soildata = [np.array(indata.MU_GLOBAL),
                np.array(indata.T_USDA_TEX_CLASS,dtype=np.int32),
                np.array(indata.T_BULK_DENSITY,dtype=np.float),
                np.array(indata.T_OC,dtype=np.float),
                np.array(indata.DRAINAGE,dtype=np.int32)]
    subsoil = [np.array(indata.MU_GLOBAL),
               np.array(indata.S_USDA_TEX_CLASS,dtype=np.int32),
               np.array(indata.S_BULK_DENSITY,dtype=np.float),
               np.array(indata.S_OC,dtype=np.float)]
    
    infiles = [os.path.join(__location__,basinMask),
               os.path.join(__location__,HWSD),
               os.path.join(__location__,basinElv),
               os.path.join(__location__,AnnPrecip),
               os.path.join(__location__,Slope)]
    
    for i in range(len(infiles)):
        
        ds = gdal.Open(infiles[i],GA_ReadOnly)
        b1 = ds.GetRasterBand(band)
        var = BandReadAsArray(b1)
        
        if i == 0:
            data = np.zeros((ds.RasterYSize,ds.RasterXSize, len(infiles)))
            gt = ds.GetGeoTransform()
            lon0 = gt[0] + (gt[1] / 2.)
            lon1 = gt[0] + (data.shape[1]*gt[1])
            lat0 = gt[3] + (data.shape[0]*gt[-1])
            lat1 = gt[3] + (gt[-1] /2)
            
        if i == 2:
            NoData = b1.GetNoDataValue()
        #
        data[:,:,i] = var[:,:]
            
        ds = None
        b1 = None
        
    lons = np.linspace(lon0,lon1,data.shape[1])
    lats = np.linspace(lat0,lat1,data.shape[0])
    
    xx,yy = np.meshgrid(lons,lats)
    
    yy = np.flipud(yy)
    
    soilfile = os.path.join(__location__,outsoil)

    if os.path.exists(soilfile)==True:
        os.remove(soilfile)
    
    cells = 0
    
    with open(soilfile, 'a') as f:
        
        cnt = 1
    
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                
                run = int(data[i,j,0])
                if run <= 0:
                    run=0       
                if data[i,j,2] == NoData:
                    run = 0
                
                if run == 0:
                    cnt+=1
                    

                else:
                    hwsdcls = data[i,j,1]
                    #print scls
                    soildic = get_usda_cls(hwsdcls, soildata,subsoil)
                    cells+=1
                               
                    soilDrain = drainAttributes[soildic['drainage']-1]['properties']
                    topSoilPro = soilAttributes[soildic['topUSDA']-1]['properties']
                    subSoilPro = soilAttributes[soildic['subUSDA']-1]['properties']
                    
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
                    
                    grdc = cnt
                    lat = yy[i,j]
                    lon = xx[i,j]
                    infilt = b_val
                    Ds = Ds_val
                    Dsmax = (data[i,j,4]/100.) * (float(subSoilPro['SatHydraulicCapacity'])*240)
                    Ws = Ws_val
                    c = 2
                    expt = 3+(2*float(topSoilPro['SlopeRCurve']))
                    expt1 = 3+(2*float(subSoilPro['SlopeRCurve']))
                    tksat = (float(topSoilPro['SatHydraulicCapacity'])*240)
                    sksat = (float(subSoilPro['SatHydraulicCapacity'])*240)
                    phis = -999
                    elev = data[i,j,2]
                    depth = 0.10
                    depth1 =s2
                    depth2 =s3
                    avg_t = 27
                    dp = 4
                    tbub = topSoilPro['BubblingPressure']
                    sbub = subSoilPro['BubblingPressure']
                    quartz = topSoilPro['Quartz']
                    quartz1 = subSoilPro['Quartz']
                    bulk_den = float(soildic['topBulkDen'])
                    bulk_den1 =float(soildic['subBulkDen'])
                    soil_den = 2650.
                    soil_den1 = 2685.
                    t_oc = float(soildic['topOC'])
                    s_oc = float(soildic['subOC'])
                    org_bulk_den = 0.25*bulk_den
                    org_bulk_den1 = 0.25*bulk_den1
                    org_soil_den = 1295.
                    org_soil_den1 = 1300.
                    off_gmt = lon * 24 / 360.
                    wrc_frac = (float(topSoilPro['FieldCapacity'])/
                                float(topSoilPro['Porosity']))
                    wrc_frac1 = (float(subSoilPro['FieldCapacity'])/
                                float(subSoilPro['Porosity']))
                    wpwp_frac = (float(topSoilPro['WiltingPoint'])/
                                float(topSoilPro['Porosity']))
                    wpwp_frac1 = (float(subSoilPro['WiltingPoint'])/
                                float(subSoilPro['Porosity']))
                    rough = 0.01
                    srough = 0.001
                    annprecip = data[i,j,3]
                    resid = topSoilPro['Residual']
                    resid1 = subSoilPro['Residual']
                    fs_act = 1
                    init_moist = (bulk_den / soil_den) * depth *1000
                    initmoist2 = (bulk_den1 / soil_den1) * depth1 *1000
                    initmoist3 = (bulk_den1 / soil_den1) * depth2 *1000
                    
                    f.write('{0}\t{1}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8}\t{9}\t{10}\t{10}\t{11}\t{12}\t{12}\t{13}\t{13}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{24}\t{25}\t{26}\t{26}\t{27}\t{28}\t{28}\t{29}\t{30}\t{30}\t{31}\t{32}\t{32}\t{33}\t{34}\t{34}\t{35}\t{36}\t{36}\t{37}\t{38}\t{39}\t{39}\t{40}\t{41}\t{41}\t{42}\t{43}\t{44}\t{45}\t{46}\t{46}\t{47}\n'.format(run,
                                        grdc,lat,lon,infilt,Ds,Dsmax,Ws,c,expt,expt1,tksat,sksat,
                                        phis,init_moist,initmoist2,initmoist3,elev,depth,
                                        depth1,depth2,avg_t,dp,tbub,sbub,quartz,quartz1,bulk_den,
                                        bulk_den1,soil_den,soil_den1,t_oc,s_oc,org_bulk_den,org_bulk_den1,
                                        org_soil_den,org_soil_den1,off_gmt,wrc_frac,wrc_frac1,wpwp_frac,
                                        wpwp_frac1,rough,srough,annprecip,resid,resid1,fs_act))
                                                                        
                    cnt+=1
                
    return

def main():
    format_soil_params(sys.argv[1],sys.argv[2],sys.argv[3],
                       sys.argv[4],sys.argv[5],sys.argv[6])

    return


if __name__ == "__main__":
    main()
