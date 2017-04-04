#------------------------------------------------------------------------------
# FILE: format_flux_climo.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 28 Oct. 2016
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script takes grided daily flux data from the VIC and
#          calculates the monthly and yearly climatological and accummulation
# DEPENDENCIES: numpy, pandas, netCDF4, osgeo (gdal)
#------------------------------------------------------------------------------

# import dependencies
import os
import sys
import netCDF4
import numpy as np
import datetime

def format_flux_climo(influx):
    
    # set up lists and dictionaries to look up month information for leap and non-leap years
    mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    nidx = {'Jan':['01',0,31],'Feb':['02',31,59],'Mar':['03',59,90],'Apr':['04',90,120],
        'May':['05',120,151],'Jun':['06',151,181],'Jul':['07',181,212],'Aug':['08',212,243],
        'Sep':['09',243,273],'Oct':['10',273,304],'Nov':['11',304,334],'Dec':['12',334,364]}
        
    lidx = {'Jan':['01',0,31],'Feb':['02',31,60],'Mar':['03',60,91],'Apr':['04',91,121],
        'May':['05',121,152],'Jun':['06',152,182],'Jul':['07',182,213],'Aug':['08',213,244],
        'Sep':['09',244,274],'Oct':['10',274,305],'Nov':['11',305,335],'Dec':['12',335,365]}
        
    
    # get the file name
    filet1 = influx
    
    # open flux file for reading
    nc = netCDF4.Dataset(filet1)
    draw = nc.variables[filet1.split('/')[-1].split('_')[0]]
    var = filet1.split('/')[-1].split('_')[0] # get variable name
    
    # get lat/lon arrays
    clon = nc.variables['longitude'][:]
    clat = nc.variables['latitude'][:]
    
    # find timing information
    toff = nc.variables['time'].units.split(' ')[-1].split('-')
    t1 = datetime.date(int(toff[0]),int(toff[1]),int(toff[2]))
    t2 = t1 + datetime.timedelta(np.max(nc.variables['time']))
            
    # get number of years in series
    yrs = np.arange(t1.year,t2.year+1)
        
    dt = -365 # day counter
    
    # set blank arrays to save data to
    outdata = np.zeros([yrs.size,len(mons),clat.size,clon.size])
    outmon = np.zeros([len(mons),clat.size,clon.size])
    
    yrvars = np.zeros([yrs.size,clat.size,clon.size]) # ???
        
    # loop over all of the years
    for y in range(yrs.size):
        # get number of days for leap and non-leap years
        if yrs[y]%4 == 0:
            idx = lidx
            dt+=366
        else:
            idx = nidx
            dt+=365
            
        # loop over each month        
        for m in range(len(mons)):
            
            midx = idx[mons[m]]
        
            t1 = midx[1] + dt
            t2 = midx[2] + dt
            
            # calculate the accumulated flux
            outdata[y,m,:,:] = np.sum(draw[t1:t2,:,:],axis=0)
            
    # mask the no data values
    outdata = np.ma.masked_where(outdata<=0,outdata)
            
    # do some accumlations over axises
    outyr = np.sum(outdata,axis=1) # yearly accumulation
    
    outmon = np.mean(outdata,axis=0) # monthly climatology
    
    avgyr = np.mean(outyr,axis=0) # yearly climatological average
    
    # get path to outfile 
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__)))
    
    # output file name
    filepath = r'{0}{1}_climatology.nc'.format(filet1.split('{0}'.format(var))[0],var)
    
    # join paths
    outfile = os.path.join(__location__,filepath)

    # open output netCDF file for writing    
    ncfile = netCDF4.Dataset(outfile, "w")
    
    # write netCDF metadata
    ncfile.Conventions = "CF-1.6"
    ncfile.title = "VIC hydrologic flux climatology"
    ncfile.source = 'VIC hydrologic model 4.2.d'
    ncfile.history = "Created using the NASA SERVIR VICUtils package. " + datetime.date.today().isoformat()
    ncfile.date_created = str(datetime.datetime.now())
    ncfile.references = "N/A"
    ncfile.comment = "N/A"
    
    #create dimensions
    ncfile.createDimension("longitude", clon[:].size)
    ncfile.createDimension("latitude", clat[:].size)
    ncfile.createDimension("Year",yrs.size)
    ncfile.createDimension("Month",len(mons))
    
    #create variables
    latvar = ncfile.createVariable("latitude", float, ("latitude",))
    latvar.long_name = "Latitude"
    latvar.units = "degrees_north"
    latvar[:] = clat[:]
    
    lonvar = ncfile.createVariable("longitude", float, ("longitude",))
    lonvar.long_name = "Longitude"
    lonvar.units = "degrees_east"
    lonvar[:] = clon[:]
    
    timevar = ncfile.createVariable("mon", float, ("Month",))
    timevar.long_name = "Month"
    timevar.units = "Calendar Month"
    timevar[:] = range(len(mons))
    
    yrvar = ncfile.createVariable("yr", float, ("Year",))
    yrvar.long_name = "Year"
    yrvar.units = "Year since 0000"
    yrvar[:] = yrs[:]
    
    data_rvar = ncfile.createVariable('monraw', float, ("Year","Month","latitude","longitude"))
    data_rvar.long_name = 'Monthly accumulated {0} for each year'.format(var)
    data_rvar.missing_value = -9999.
    data_rvar.units = "mm/mon"
    data_rvar[:] = outdata[:,:,:,:]
    
    data_yvar = ncfile.createVariable('yrraw', float, ("Year","latitude","longitude"))
    data_yvar.long_name = 'Yearly accumulated {0} for each year'.format(var)
    data_yvar.missing_value = -9999.
    data_yvar.units = "mm/yr"
    data_yvar[:] = outyr[:,:,:]
    
    data_mvar = ncfile.createVariable('monclimo', float, ("Month","latitude","longitude"))
    data_mvar.long_name = 'Monthly accumulated {0} climatology'.format(var)
    data_mvar.missing_value = -9999.
    data_mvar.units = "mm/mon"
    data_mvar[:] = outmon[:,:,:]
    
    data_avar = ncfile.createVariable('yrclimo', float, ("latitude","longitude"))
    data_avar.long_name = 'Yearly accumulated {0} climatology'.format(var)
    data_avar.missing_value = -9999.
    data_avar.units = "mm/yr"
    data_avar[:] = avgyr[:,:]
    
    # close netCDF file
    ncfile.close()
    
    print 'Flux climotology file written to "{0}" directory'.format(__location__)

    return 
    
def main():
    # checking user input
    if len(sys.argv) != 2:
        print "Wrong user input"
        print "Caluculate flux climatology from flux NetCDF"
        print "usage: python format_flux_climo.py <vic flux netcdf file>"
        #print "DIR INPUTS SHOULD CONTAIN TRAILING /"
        sys.exit()
        
    else:
        # pass system arguments to the function
        format_flux_climo(sys.argv[1])
        
    return
    
# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()