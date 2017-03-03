import netCDF4
import numpy as np
import datetime

def find_nearest_idx(arr,var):    
    return (np.abs(arr-var)).argmin()
    
bb = [0.125,34.75,-0.5,36.0]
    
inpath = '/Users/kmarkert/Documents/Kel/UAH/classes/ESS/thesis/models/input/forcing/chirps_climo/chirps-v2.0.YYYY.days_p05.nc'

yrs = np.arange(2013,2014)

for i in range(len(yrs)):
    infile = inpath.replace('YYYY',str(yrs[i]))
    
    nc = netCDF4.Dataset(infile)
    
    if i == 0:
        
        lats = nc.variables['latitude'][:]
        lons = nc.variables['longitude'][:]
                
        img_bb = [find_nearest_idx(lats,bb[0]),find_nearest_idx(lons,bb[1]),
                  find_nearest_idx(lats,bb[2]),find_nearest_idx(lons,bb[3])]
        print img_bb
                
        toff = datetime.date(1980,1,1)
        inidate = toff + datetime.timedelta(int(np.min(nc.variables['time'])))
        enddate = toff + datetime.timedelta(int(np.max(nc.variables['time'])))
        
        if yrs[i]%4 == 0:
            dayOff = 2
        else:
            dayOff = 1
        
        days = (enddate-inidate).days +1
        
    else:
        pass
        
    pr = nc.variables['precip'][:,:,:]
    
    nc.close()
    
    outpr = pr[:,img_bb[2]:img_bb[0],img_bb[1]:img_bb[3]]
    outlats = lats[img_bb[2]:img_bb[0]]
    outlons = lons[img_bb[1]:img_bb[3]]
    
    pr = None
    
    ncfile = netCDF4.Dataset(infile[:-3]+"_NyandoBasin.nc", "w")

    ncfile.Conventions = "CF-1.6"
    ncfile.title = "CHIRPS Version 2.0"
    ncfile.source = 'http://chg.geog.ucsb.edu/data/chirps/index.html'
    ncfile.history = "created by Climate Hazards Group, clipped for smaller scale analysis. " + datetime.date.today().isoformat()
    ncfile.date_created = str(datetime.datetime.now())
    ncfile.references = "Funk, C.C., Peterson, P.J., Landsfeld, M.F., Pedreros, D.H., Verdin, J.P., Rowland, J.D., Romero, B.E., Husak, G.J., Michaelsen, J.C., and Verdin, A.P., 2014, A quasi-global precipitation time series for drought monitoring: U.S. Geological Survey Data Series 832, 4 p., http://dx.doi.org/110.3133/ds832. "
    ncfile.comment = "N/A"
    
    ncfile.start_date = inidate.isoformat()
    ncfile.end_date = enddate.isoformat()
    
    #create dimensions
    ncfile.createDimension("longitude", len(outlons))
    ncfile.createDimension("latitude", len(outlats))
    ncfile.createDimension("time", days)
    
    #create variables
    latvar = ncfile.createVariable("latitude", float, ("latitude",))
    latvar.long_name = "Latitude"
    latvar.units = "degrees_north"
    latvar[:] = outlats
    
    lonvar = ncfile.createVariable("longitude", float, ("longitude",))
    lonvar.long_name = "Longitude"
    lonvar.units = "degrees_east"
    lonvar[:] = outlons
    
    timevar = ncfile.createVariable("time", int, ("time",))
    timevar.long_name = "Time"
    timevar.units = "days since " + inidate.isoformat()
    timevar.calendar = 'gregorian'
    timevar[:] = range(0, days)
    
    data_var = ncfile.createVariable('precip', float, ("time","latitude","longitude"))
    data_var.long_name = "Climate Hazards group InfraRed Precipitation with Stations"
    data_var.missing_value = -9999.0
    data_var.units = "mm"
    data_var[:] = outpr[:,:,:]
    
    outpr = None
    
    ncfile.close()