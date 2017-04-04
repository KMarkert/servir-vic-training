#!/usr/bin/python

#----------------------------------------------------
# Program to convert VIC fluxes files to NetCDF file
# will ask the user wich variable he wants to export
# and also for wich years. Assumes there is data
# for the entire time period, from 1-jan to 31-dec
# SET UP FOR DAILY TIME STEP. FLUX FILE SHOUD NOT
# CONTAIN HOUR RECORD!!
#----------------------------------------------------

#------------------------------------------------
# Writen by Daniel de Castro Victoria
# dvictori@cena.usp.br or daniel.victoria@gmail.com
# Needs python libraries Numeric and Scientific
# 03-dec-2004
#
# Script updated by Kel Markert
# kel.markert@nasa.gov or kel.markert@uah.edu
#-------------------------------------------------

# import dependencies
import sys
import os, string
# handle dates...
import datetime as dt
# NetCDF and Numeric
from netCDF4 import *
from numpy import *

def flux2nc(influxes,outpath,var=None,start_year=None,end_year=None):
    
    # building file list and sorted lat lon list
    dirin = os.path.dirname(influxes)
    
    try:
        file_list = os.listdir(dirin)
    except OSError:
        raise OSError('Input flux directory not valid, please fix path')
    
    lat_t = []
    lon_t = []
    lat = []
    lon = []
    
    try:
        for f in file_list:
            lat_t.append(float(string.split(f, "_")[1]))
            lon_t.append(float(string.split(f, "_")[2]))
    except ValueError:
        raise ValueError('Input path contains files that are not flux files')
    
    for i in lat_t:
        if i not in lat:
            lat.append(i)
    
    for i in lon_t:
        if i not in lon:
            lon.append(i)
    
    
    # putting in order. Lat should be from top to botom
    # lon from left to right
    lon.sort()
    lat.sort()
    lat.reverse()
    
    del(lat_t)
    del(lon_t)
    
    # if variable is not set, get it from user
    if var == None:
    
        #determining the parameter to use
        print "Choose output parameter"
        print "1 - Precipitation"
        print "2 - Evapotranspiration"
        print "3 - Runoff"
        print "4 - Base flow"
        print "5 - Snow Water Equivalent"
        print "6 - Soil moisture"
        varini = input('Choose output (1 a 6)>')
        
        #getting the collumn right
        if varini < 6:
            var = varini + 2
        elif varini == 6:        #more than one soil layer...
            camada = input('which soil layer?>')
            var = varini + 2 + camada
    
    #set name of out_file. Named after parameter choice
    if var == 3:
        var_txt = "ppt"
        var_name = "Precipitation"
    elif var == 4:
        var_txt = "evap"
        var_name = "Evapotranspiration"
    elif var == 5:
        var_txt = "runoff"
        var_name = "Runoff"
    elif var == 6:
        var_txt = "base"
        var_name = "Baseflow"
    elif var == 7:
        var_txt = "swe"
        var_name = "Snow Water Equivalent"
    else:
        var_txt = "soilLyr"+str(camada)
        var_name = "Soil moisture, layer {0}".format(camada)
    
    # if the date information is not set get it from user
    if start_year == None:
        # for what date?
        start_year = input("Enter start year:")
    if end_year == None:
        end_year = input("End year:")
    
    # set date information in datetime object
    inidate = dt.date(start_year,01,01)
    enddate = dt.date(end_year,12,31)
    
    # calculate number of days in time series
    days = enddate.toordinal() - inidate.toordinal()+1
    
    #print "Gridding {0} data...".format(var_name)
    
    #
    # create array containig all data
    # This is going to be huge. Create an array with -9999 (NoData)
    # Then populate the array by reading each flux file
    #
    
    all_data = zeros([days,len(lat),len(lon)], dtype=float32)
    all_data[:,:,:] = -9999
    
    c = len(file_list)
    
    # for each file in list
    for f in file_list:
        # get lat & lon and it's index
        latitude = float(string.split(f, sep="_")[1])
        longitude = float(string.split(f, sep="_")[2])
        lat_id = lat.index(latitude)
        lon_id = lon.index(longitude)
    
        c = c -1
        
        infile = open(dirin+'/'+f, "r")
        lixo = infile.readlines()
        infile.close()
        dado = []
    
        for l in lixo:
            if int(string.split(l, sep="\t")[0]) in range(inidate.year, enddate.year+1):
                dado.append(float(string.split(l, sep="\t")[var]))
            # putting data inside array.
            # Since data has lat & lon fixed uses dimension [:,lat_index,lon_index]
    
        all_data[:,lat_id,lon_id] = dado
    
    del dado # del data to free memory for large datasets
    
    try:
    
        # open netCDF file for writing
        ncfile = Dataset(outpath+str(var_txt)+'_'+str(start_year)+".nc", "w")
        
        # set netCDF metadata information
        ncfile.Conventions = "CF-1.6"
        ncfile.title = "VIC hydrologic flux outputs"
        ncfile.source = 'VIC hydrologic model 4.2.d'
        ncfile.history = "Created using the script created by NASA SERVIR. " + dt.date.today().isoformat()
        ncfile.date_created = str(dt.datetime.now())
        ncfile.references = "N/A"
        ncfile.comment = "N/A"
        
        ncfile.start_date = inidate.isoformat()
        ncfile.end_date = enddate.isoformat()
        
        #create dimensions
        ncfile.createDimension("longitude", len(lon))
        ncfile.createDimension("latitude", len(lat))
        ncfile.createDimension("time", days)
        
        #create variables
        latvar = ncfile.createVariable("latitude", float, ("latitude",))
        latvar.long_name = "Latitude"
        latvar.units = "degrees_north"
        latvar[:] = lat
        
        lonvar = ncfile.createVariable("longitude", float, ("longitude",))
        lonvar.long_name = "Longitude"
        lonvar.units = "degrees_east"
        lonvar[:] = lon
        
        timevar = ncfile.createVariable("time", int, ("time",))
        timevar.long_name = "Time"
        timevar.units = "days since " + inidate.isoformat()
        timevar.calendar = 'gregorian'
        timevar[:] = range(0, days)
        
        # save gridded flux data to file
        data_var = ncfile.createVariable(var_txt, float, ("time","latitude","longitude"))
        data_var.long_name = var_name
        data_var.missing_value = -9999.0
        data_var.units = "mm"
        data_var[:] = all_data[:,:,:]
        
        # close the file
        ncfile.close()
        
    except IOError:
        raise IOError('Output path is not valid, please fix the path string')
    
    return

def main():
    # checking user input
    if len(sys.argv) != 3:
        print "Wrong user input"
        print "Convert VIC fluxes files to NetCDF"
        print "usage flux2cdf.py <vic flux dir> <out netcdf dir>"
        print "DIR INPUTS SHOULD CONTAIN TRAILING /"
        sys.exit()
    
    if sys.argv[1][-1] != "/":
        print "VIC FLUX DIR SHOULD CONTAIN TRAILING /"
        print "fixing it for you..."
        sys.argv[1] = sys.argv[1] + "/"
    
    print "IMPORTANT: "+sys.argv[1]+" SHOULD CONTAIN ONLY FLUXES FILES!!!"
    
    flux2nc(sys.argv[1],sys.argv[2])
    
    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    main()
