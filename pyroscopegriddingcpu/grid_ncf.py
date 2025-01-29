# Saving as NetCDF
# 
# There are several ways to save gridded data as netcdf:
# 
# 1. Translation from geoTIFF format to netcdf 
# - this is taken care of by gdal
# - does not support time dimension
# 2. Single file saving: gridNC_single(limit, gsize, indata, inlat, inlon, filename)
# - limit - [lat1, lat2, lon1, lon2]
# - gsize - pixel size
# - indata - raw satellite data from import in
# - inlat
# - inlon
# - filename (that we save to, including path)
# 3. Single time file saving: gridNC_time(limit, gsize, indata, inlat, inlon, filename)
# - limit - [lat1, lat2, lon1, lon2]
# - gsize - pixel size
# - indata - raw satellite data from import in (multiple file array)
# - inlat
# - inlon
# - filename (that we save to, including path)
# 4. Multi time file saving: 

from http.client import MOVED_PERMANENTLY
import netCDF4
import numpy as np
import numpy.ma as ma
import time
import os

"""
from gridding import *
from sat_data_input import *
from filter_data import *
from naming_conventions import *
from solar_zenith import *
from fileparser import *
"""
from pyroscopegriddingcpu.gridding import *
from pyroscopegriddingcpu.sat_data_input import *
from pyroscopegriddingcpu.filter_data import *
from pyroscopegriddingcpu.naming_conventions import *
from pyroscopegriddingcpu.solar_zenith import *
from pyroscopegriddingcpu.fileparser import *

#import gridding
#import sat_data_input
#import filter_data
#import naming_conventions
#import solar_zenith

import datetime

# full satellite list 
# used to check if any satellites are missing
#full_satellite_list = ['AERDT_L2_ABI_G16', 'AERDT_L2_ABI_G17', 'AERDT_L2_AHI_H08', 'AERDT_L2_VIIRS_SNPP', 'MOD04_L2', 'MYD04_L2']
#full_satellite_list = ['ABI_G16', 'ABI_G17', 'AHI_H08', 'AHI_H09','VIIRS_NOAA20','VIIRS_SNPP', 'MOD04', 'MYD04']
#full_satellite_list = ['ABI_G16', 'ABI_G17', 'AHI_H08', 'VIIRS_SNPP', 'MOD04', 'MYD04']
full_satellite_list = ['MOD04', 'MYD04', 'VIIRS_SNPP', 'VIIRS_NOAA20', 'ABI_G16', 'ABI_G17', 'AHI_H08', 'AHI_H09']

# takes in dictionary 
# keys are satellite names, values = array of filepaths
# returns 1D array of satellite file inputs
def sat_list_name_concat(file_list):
    file_inputs = []
    
    for key in file_list:
        for f in file_list[key]:
            f = convert_path_to_linux(f)
            file_inputs.append(f[(f.rfind("/") + 1):])
    
    return file_inputs


# filename - 'gridNet.nc' (this is what it saves to)
# sat_files - list of satellite files to pull data from
# saves all satellite sensor gridded data as separate variables
#def grid_nc_single_statistics(limit, gsize, inlat, inlon, geo_list, phy_list, filelist, filename):
def grid_nc_sensor_statistics_metadata(limit, gsize, geo_list, phy_list, filelist, filename, time_start, time_diff, 
                                       phy_nc=None, phy_hdf=None, static_file = None, pixel_range = None):
    #start_timer = time.time()
    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    lat_dim = int(1+round(((maxlat-minlat)/gsize))) #int(1+round(((maxlat-minlat)/gsize)))
    lon_dim = int(1+round(((maxlon-minlon)/gsize))) #int(1+round(((maxlon-minlon)/gsize)))
    
    #set target netcdf file
    ds = netCDF4.Dataset(filename, 'w',format='NETCDF4')
    
    #start_timer = time.time()
    time_date = time_start.strftime('%Y-%m-%d')
    time_time = time_start.strftime('%H:%M:%S.%f')[:-4]
    time_end = time_start + datetime.timedelta(minutes=time_diff)
    time_end_date = time_end.strftime('%Y-%m-%d')
    time_end_time = time_end.strftime('%H:%M:%S.%f')[:-4]

    #set global attributes
    ds.title = "Level-3 quarter-degree 30-minute global aerosol data gridded, averaged and merged from LEO and GEO sensors"
    ds.references = "1) Levy, R. C., S. Mattoo, L. A. Munchak, et al. 2013. The Collection 6 MODIS Aerosol Products over Land and Ocean. Atmos Meas Tech 6 2989-3034 [10.5194/amt-6-2989-2013]; 2) Gupta, P.; Remer, L.A.; Patadia, F.; Levy, R.C.; Christopher, S.A. High-Resolution Gridded Level 3 Aerosol Optical Depth Data from MODIS. Remote Sens. 2020, 12, 2847. https://doi.org/10.3390/rs12172847"
    ds.history = "pyroscopegriddingcpu version 1.5.0.0" 
    ds.institution = "NASA Goddard Space Flight Center, Climate and Radiation Laboratory"
    ds.publisher_name = "LAADS"
    ds.publisher_url = "https://ladsweb.modaps.eosdis.nasa.gov"
    ds.LongName = "THANGS Combined GEOLEO Aerosol 30-Min average L3 Global 0.25x0.25 degree grid"
    ds.ShortName = "XAERDT_L3_MEASURES_QD_HH"
    ds.VersionID = "001"
    ds.identifier_product_doi = "10.5067/MEaSUREs/GLDTA/XAERDT_L3_MEASURES_QD_HH.001"
    ds.identifier_product_doi_authority  = "https://dx.doi.org/"
    ds.Conventions = "CF-1.8"
    ds.Format = "NetCDF-4"
    ds.processing_level = "L3"
    ds.AlgorithmType = "SCI"
    ds.platform = "Terra/Himawari-8,9/Aqua//NOAA-20/GOES-16,17/SNPP"
    ds.instrument = "MODIS+AHI+ABI+VIIRS"
    ds.SatelliteInstrument = "THANGS"

    filename = convert_path_to_linux(filename) #convert to linux pwd
    cut_index = filename.rindex("/")
    filename_without_path = filename[cut_index+1 :]
    ds.LocalGranuleID = filename_without_path
    
    #production_yyyymmdd = filename_without_path[-11:-3]
    production_yyyymmdd = '{}-{}-{}'.format(filename_without_path[-11:-7],
            filename_without_path[-7:-5],filename_without_path[-5:-3])
    currtime = datetime.datetime.now()
    production_hhmmss = currtime.strftime("%H:%M:%S")
    ds.RangeBeginningDate = time_date
    ds.RangeBeginningTime = time_time
    ds.RangeEndingDate = time_end_date
    ds.RangeEndingTime = time_end_time
    ds.ProductionTime = production_yyyymmdd +" "+ production_hhmmss
    ds.geospatial_lon_units = "degrees_east" 
    ds.WestBoundingCoordinate = "-179.875" 
    ds.EastBoundingCoordinate = "179.875" 
    ds.geospatial_lat_units = "degrees_north" 
    ds.NorthBoundingCoordinate = "89.875" 
    ds.SouthBoundingCoordinate = "-89.875" 
    
    
    #list of all files used in creation of this output
    ds.InputPointer = ','.join(sat_list_name_concat(filelist))
    ds.ancillary_files = "LSM_ELV_QDEG_FIXED.nc" 
    ds.PGE_StartTime = time_date + " " + time_time  
    ds.PGE_EndTime =  time_end_date + " " + time_end_time
    ds.PGE_Name = "PGE105"
    ds.PGEVersion = "1.0.0-1" 
    
    # Create time, lat, lon, sensor dimensions
    timed = ds.createDimension('time', None)
    lat = ds.createDimension('lat', lat_dim) #latitude is fatitude (y)
    lon = ds.createDimension('lon', lon_dim) # x
    #sensor = ds.createDimension('sensor', 6) #only used for SensorWeighting mask
    sensor = ds.createDimension('sensor', len(full_satellite_list))
    
    #set time variable to filename parse
    times = ds.createVariable('time', 'f4', ('time',)) #
    times.long_name = "time"
    times.calendar = "standard"
    times.units = "minutes since "+str(time_start)
    
    times[:] = [0.]

    
    #set lats and lons  (or rather y and x)
    lats = ds.createVariable('lat', 'f4', ('lat',), fill_value=-9999.,compression='zlib') #latitude is fatitude
    lons = ds.createVariable('lon', 'f4', ('lon',), fill_value=-9999.,compression='zlib')
    sensor_var = ds.createVariable("sensors", np.short, ('sensor', ))
    
    #sensors metadata
    sensor_var.long_name = "Sensors: 1- MODIS-Terra, 2 - MODIS-Aqua, 3 - VIIRS-SNPP, 4 - VIIRS-NOAA20, 5 - ABI-G16, 6 - ABI-G17, 7 - AHI-H08, 8 - AHI-H09"
    #sensor_var.units = "None"
    
    #latitude / longitude metadata
    lats.valid_range = [-90.0, 90.0]
    lats.standard_name = "latitude"
    lats.long_name = "Geodectic Latitude"
    lats.units = "degrees_north"
    lats.scale_factor = 1.
    lats.add_offset = 0.
    lats.Parameter_Type = "Equal angle grid center location"
    lats.Geolocation_Pointer = "Geolocation data not applicable" 
    lats._CoordinateAxisType = "Lat"
    
    lons.valid_range = [-180.0, 180.0]
    lons.standard_name = "longitude"
    lons.long_name = "Geodectic Longitude"
    lons.units = "degrees_east"
    lons.scale_factor = 1.
    lons.add_offset = 0.
    lons.Parameter_Type = "Equal angle grid center location"
    lons.Geolocation_Pointer = "Geolocation data not applicable" 
    lons._CoordinateAxisType = "Lon"

    #bounds for lats and lons
    lats[:] = np.arange(0,lat_dim)*gsize+minlat
    lons[:] = np.arange(0,lon_dim)*gsize+minlon
    
    # unknown number of sensors, save values to list
    values = []
    sensors = {}

    # instantiate empty arrays for summary statistics
    meta = []
    indata_temp = []
    inlat_temp = []
    inlon_temp = []
    
    # copy static file data into this file
    # Land_Sea_Mask, Topographic_Altitude
    if(static_file != None):
        static_vars = ["Land_Sea_Mask", "Topographic_Altitude"]
        
        #with netCDF4.Dataset(static_file) as src, ds as dst:
        src = netCDF4.Dataset(static_file)
        # copy dimensions
        """
        for name, dimension in src.dimensions.items():
            ds.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        """
        # copy all file data in the specified variables
        for name, variable in src.variables.items():
            if name in static_vars:
                #kvp
                #dimensions = ('lat', 'lon',)
                dimensions = ('time', 'lat', 'lon',)
                x = ds.createVariable(name, variable.datatype, dimensions,compression='zlib') #variable.dimensions)
                #ds[name][:, :] = src[name][:]
                ds[name][0, :, :] = src[name][:]
                # copy variable attributes all at once via dictionary
                ds[name].setncatts(src[name].__dict__)
                
                # add manual metadata
                if name == "Topographic_Altitude":
                    ds[name].long_name = "Averaged topographic altitude (in meter) for Land"
                    ds[name].units = "m"
                    #kvp
                    ds[name].valid_range = [-1, 10000]
                    ds[name].scale_factor = 1.0
                    ds[name].add_offset = float(0)

                else:
                    #kvp
                    #ds[name].long_name = "Land_Sea_Flag (based on MOD03 Landsea mask 0 = Ocean, 1 = Land and ephemeral water 2 = Coastal)"
                    ds[name].long_name = "Land_Sea_Flag (Landsea mask: 0 = Ocean, 1 = Land and ephemeral water)"
                    ds[name].valid_range = [0, 1]
                    #ds[name].scale_factor = 1.0
                    #ds[name].add_offset = 0.0

        src.close()
    
    print("Done with static file transfer")
    
    # Solar zenith angle 
    name = "Solar_Zenith_Angle"
    #kvp
    #solar_zenith_variable = ds.createVariable(name, 'f4', ('lat', 'lon', ), fill_value=int(-9999),compression='zlib')
    solar_zenith_variable = ds.createVariable(name, 'f4', ('time', 'lat', 'lon', ), fill_value=int(-9999),compression='zlib')
    solar_zenith_variable.long_name = "Solar Zenith Angle at the center of grid calculated for central time of the file"
    solar_zenith_variable.units = "degrees"
    solar_zenith_variable.valid_range = [ 0, 18000]
    solar_zenith_variable.scale_factor = 0.01
    solar_zenith_variable.add_offset = float(0)
    #netCDF4.nc_put_att(ds, solar_zenith_variable, 'add_offset', 'f4', 1, 0)
    #kvp
    #solar_zenith_variable[:, :] = get_SZA_parallelized(limit, gsize, time_start, time_diff)
    solar_zenith_variable[0, :, :] = get_SZA_parallelized(limit, gsize, time_start, time_diff)
    print("solar zenith done")

    # for calculating LEOGEO statistics
    # LEOGEO: Mean, STD, NumberOfSensors, SensorWeighting, TotalPixels
    # index of value corresponds to individual sensor
    # leogeo_stats = {"Mean": [], "STD": [], "TotalPixels":[]}
    leogeo_stats_arr = {} #contains leogeo_stats as value, geophys as key
    # sensor_idx_stats = {"NumberOfSensors":[], "SensorWeighting":[]}
    sensor_idx_arr = {} #key = geophys value, value = sensor_idx
    
    # instantiate the satellites that are not present
    # no inputs for these satellites, blank data for their netCDF4 output variables
    # not_present_satellites = [s for s in full_satellite_list if s not in list(filelist.keys())]
    # accounts for abbreviations
    present_satellites = []
    for s in full_satellite_list: #check all satellites
        for present_s in list(filelist.keys()):
            if s in present_s: #satellite is present
                present_satellites.append(s)
                
    not_present_satellites = [s for s in full_satellite_list if s not in present_satellites]
    #kvp
    print("Present Sensors: ",present_satellites) 
    print('Sensors not available: ',not_present_satellites) 
    # metadata for not_present_satellites

    for s_name in not_present_satellites: # not present satellites
        for j, p_vars in enumerate(phy_list): # creating a variable for each geophysical variable 
            if not("Sensor_Zenith_Angle" in p_vars or "Scattering_Angle" in p_vars):
                aod_statistics = ["Mean", "STD", "Pixels"]
                
                for aod_stat in aod_statistics:
                    name = nc_var_name(p_vars, s_name, aod_stat) #str(s_name+"_"+p_vars)
                    
                    #default fill value is fill value for first satellite encountered when reading files in
                    if len(meta) > 0:
                        #kvp
                        #ds.createVariable(name, np.short, ('lat', 'lon', ), fill_value = meta[0]["_FillValue"],compression='zlib')
                        ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value = meta[0]["_FillValue"],compression='zlib')
                    else:
                        ##ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value = -9999)
                        #name=ds.createVariable(name, np.short, ('lat', 'lon', ), fill_value = -9999,compression='zlib')
                        name=ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value = -9999,compression='zlib')

                        name.units = "None"
                        for_longname = nc_long_name(p_vars, s_name, aod_stat, aod_long=None)

                        if not(aod_stat in "Pixels") or not("Pixels" in aod_stat):
                            name.valid_range = [-100, 5000]

                            if (p_vars == 'Optical_Depth_Land_And_Ocean'):
                                longname = for_longname + (" AOT at 0.55 micron for both ocean (Average) (Quality flag = 1, 2, 3) and land (corrected) (Quality flag = 3) for the grid")
                            else:
                                longname = for_longname + (" AOT at 0.55 micron for both ocean (Average) and land (corrected) with all quality data (Quality flag = 0, 1, 2, 3) for the grid")

                            name.long_name = longname
                            name.scale_factor = 0.001
                            name.add_offset = 0.0

                        else:
                            name.valid_range = [0, 10000]
                            # valid range for Pixel must be from 0
                            #if pixel_range == None:
                            if (p_vars == 'Optical_Depth_Land_And_Ocean'):
                                name.long_name =for_longname + (" AOT at 0.55 micron for both ocean (Average) (Quality flag = 1, 2, 3) and land (corrected) (Quality flag = 3) for the grid")
                            else:
                                name.long_name=for_longname + (" AOT at 0.55 micron for both ocean (Average) and land (corrected) with all quality data (Quality flag = 0, 1, 2, 3) for the grid")

                        name.Parameter_Type = "Output"

            else: #phy_var is NOT aod
                name = nc_var_name(p_vars, s_name) #str(s_name+"_"+p_vars)
                if len(meta) > 0:
                    #kvp
                    #values.append(ds.createVariable(name, np.short, ('lat', 'lon', ), fill_value = meta[0]["_FillValue"],compression='zlib'))
                    values.append(ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value = meta[0]["_FillValue"],compression='zlib'))
                else:
                    ##values.append(ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value = -9999))
                    #name=ds.createVariable(name, np.short, ('lat', 'lon', ), fill_value = -9999, compression='zlib')
                    name=ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value = -9999, compression='zlib')

                    name.units = "degrees" #meta[j]["units"]
                    name.valid_range = [0, 18000]
                    name.long_name = nc_long_name(p_vars, s_name)
                    name.scale_factor = 0.01
                    name.add_offset = 0.0
                    if ("Sensor_Zenith_Angle" in p_vars):
                        name.Parameter_Type = s_name + " Input"
                    if ("Scattering_Angle" in p_vars):
                        name.Parameter_Type = "Output"                    
   
    #run through and assign gridded data to sensor variable
    # filenames are in the form:
    # [[time0: [sat1], [sat2], ... ], ...[timek: [sat1], [sat2], ....]] 
    # we assume we receive a signle dict for a time interval though
    for i, s_name in enumerate(filelist.keys()): 
        start_timer = time.time()

        for s in filelist.get(s_name):
            print("WORK SENSOR: ", s_name)
            print("FILE ORIGINAL: ", str(s))
            if not os.path.exists(s):
                s = convert_to_other_system(s)
            
            L2FID,GeoID,PhyID = open_file(s)

            # check for hdf files 
            if (str(s).split(".")[-1] == "hdf"):
                geo_list = ['Latitude', 'Longitude']
                lat,lon,phy_vars, meta = filter_data_hdf(GeoID, PhyID, geo_list, phy_list, phy_hdf)
                L2FID.end()
                
            else: #netcdf
                geo_list = ['latitude', 'longitude']
                lat,lon,phy_vars, meta = filter_data_nc(GeoID, PhyID, geo_list, phy_list, phy_nc)
                L2FID.close()
            
            # statistics appending for this sensor
            # index for arrays = geophysical variable indices
            for count, p_var in enumerate(phy_list):
                
                #create empty array for geophysical variable
                if len(indata_temp) < count+1:
                    indata_temp.append([])
                    inlat_temp.append([])
                    inlon_temp.append([])
                
                indata_temp[count] = ma.concatenate([indata_temp[count], phy_vars[count]])
                inlat_temp[count] = ma.concatenate([inlat_temp[count], lat[count]])
                inlon_temp[count] = ma.concatenate([inlon_temp[count], lon[count]])
        # end sensor concatenations
        # indata = [[geovar data 1], [geovar data 2], ... , [geovar data n]]
        # grid parameters
        # individual sensor statistics
        for j, p_vars in enumerate(phy_list):
            
            #special metadata and statistic calculations for AOD
            # values that need LEOGEO statistic calculations
            if not("Sensor_Zenith_Angle" in p_vars or "Scattering_Angle" in p_vars):#"Optical_Depth_Land_And_Ocean" in p_vars:
                aod_statistics = ["Mean", "STD", "Pixels"]
                avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,float(gsize),indata_temp[j],inlat_temp[j],inlon_temp[j])
                mintau[np.isnan(mintau)]=meta[j]["_FillValue"]
                maxtau[np.isnan(maxtau)]=meta[j]["_FillValue"]
                avgtau[np.isnan(avgtau)]=meta[j]["_FillValue"]
                stdtau[np.isnan(stdtau)]=meta[j]["_FillValue"]
                
                #concatenate AOD data for leogeo stats later
                for aod_stat in aod_statistics:
                    
                    index = len(values) #i * len(phy_list) + j
                    aod_long = meta[j]["long_name"]
                    name = nc_var_name(p_vars, s_name, aod_stat) #str(s_name+"_"+p_vars)
                    #kvp
                    #values.append(ds.createVariable(name, np.short, ('lat', 'lon', ), fill_value=meta[j]["_FillValue"],compression='zlib'))
                    values.append(ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value=meta[j]["_FillValue"],compression='zlib'))
                    
                    # metadata
                    values[index].units = "None" #meta[j]["units"]
                    values[index].valid_range = meta[j]["valid_range"]
                    values[index].long_name = nc_long_name(p_vars, s_name, aod_stat, aod_long) + " for the grid" #meta[j]["long_name"]
                    values[index].Parameter_Type = meta[j]["Parameter_Type"]
                    
                    # check if pixels to add add_offset or scale_factor or not
                    if not(aod_stat in "Pixels") or not("Pixels" in aod_stat):
                        values[index].scale_factor = round(meta[j]["scale_factor"], 3)
                        values[index].add_offset = meta[j]["add_offset"]
                        values[index].valid_range = meta[j]["valid_range"]
                    else:
                        # valid range for Pixel must be from 0
                        if pixel_range == None:
                            values[index].valid_range = [0, 10000]
                        else:
                            values[index].valid_range = [0, 10000]

                    
                    # instantiate the dictionary for statistics
                    if not p_vars in leogeo_stats_arr:
                        leogeo_stats_arr[p_vars] = {"Mean": [], "STD": [], "TotalPixels":[]}
                        
                        
                    # make sure key is in dict 
                    # instantiation
                    if not p_vars in sensors:
                        sensors[p_vars] = []
                    
                    if not p_vars in sensor_idx_arr:
                        sensor_idx_arr[p_vars] = {}
                    
                    #statistics to be added
                    if aod_stat== "Mean":
                        rotated = list(zip(*avgtau))[::-1] # check why zipped
                        #print("AFTER GRIDDING: ", rotated)
                        leogeo_stats_arr[p_vars] ["Mean"].append(np.flipud(rotated)) #for leogeo calculation later
                        
                        #sensor idx (Sensor Weighting)
                        temp_a = np.flipud(rotated)
                        
                        temp_a[temp_a > -9999] = 1
                        temp_a[temp_a <= meta[j]["_FillValue"]] = 0 #get rid of this 2/1
                        temp_a[np.isnan(temp_a)] = 0
                        sensor_idx_arr[p_vars][get_sensor(s_name)] = temp_a
                                           
                         
                    elif aod_stat == "STD":
                        rotated = list(zip(*stdtau))[::-1]
                        
                        leogeo_stats_arr[p_vars]["STD"].append(np.flipud(rotated))
                        #print("AOD STD", rotated)
                        
                        #leogeo_stats["STD"].append(np.flipud(rotated)) #for leogeo calculation later
                    elif aod_stat == "Pixels":
                        rotated = list(zip(*count))[::-1]
                        
                        leogeo_stats_arr[p_vars]["TotalPixels"].append(np.flipud(rotated))
                        #leogeo_stats["TotalPixels"].append(np.flipud(rotated)) #for leogeo calculation later
                    else:
                        rotated = list(zip(*avgtau))[::-1]
                        
                        leogeo_stats_arr[p_vars]["Mean"].append(np.flipud(rotated))
                        #leogeo_stats["Mean"].append(np.flipud(rotated)) #for leogeo calculation later
                        print("STATISTIC ERROR")
                    
                    
                    #manually set to fill_value
                    rotated = np.array(rotated)
                    
                    final_input = np.flipud(rotated)#.astype(np.short)
                    final_input = final_input #/meta[j]["scale_factor"]
                    #kvp check here
                    #values[index][:, :] = final_input
                    #values[index][:, :] = values[index][:, :].astype("f4")
                    values[index][0, :, :] = final_input
                    values[index][0, :, :] = values[index][0, :, :].astype("f4")

                    # AOD print
                   
                    sensors[p_vars].append(np.flipud(rotated)) #add it to complete dataset for LEOGEO calculations
                    
                #time
                end_timer = time.time()
                print("Sensor: ",s_name ," time: ", end_timer - start_timer)
                
            else:  #phy_var is NOT AOD

                avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,float(gsize),indata_temp[j],inlat_temp[j],inlon_temp[j])
                mintau[np.isnan(mintau)]=meta[j]["_FillValue"]
                maxtau[np.isnan(maxtau)]=meta[j]["_FillValue"]
                avgtau[np.isnan(avgtau)]=meta[j]["_FillValue"]
                stdtau[np.isnan(stdtau)]=meta[j]["_FillValue"]

                index = len(values) #i * len(phy_list) + j
                name = nc_var_name(p_vars, s_name) #str(s_name+"_"+p_vars)
                #kvp
                #values.append(ds.createVariable(name, np.short, ('lat', 'lon', ), fill_value=meta[j]["_FillValue"],compression='zlib'))
                values.append(ds.createVariable(name, np.short, ('time', 'lat', 'lon', ), fill_value=meta[j]["_FillValue"],compression='zlib'))
                
                # metadata
                values[index].units = "degrees"#meta[j]["units"]
                values[index].valid_range = meta[j]["valid_range"]
                values[index].long_name = nc_long_name(p_vars, s_name)#meta[j]["long_name"]
                values[index].scale_factor = round(meta[j]["scale_factor"], 3)
                values[index].add_offset = meta[j]["add_offset"]
                values[index].Parameter_Type = meta[j]["Parameter_Type"]

                """
                # make sure key is in dict 
                # instantiation
                if not p_vars in sensors:
                    sensors[p_vars] = []
                """
                
                rotated = list(zip(*avgtau))[::-1]
                #print("ROTATED NON AOD")
                #print(np.flipud(rotated)) # print what is input into the variable
                
                #manually set to fill_value
                rotated = np.array(rotated)
                rotated[rotated > meta[j]["valid_range"][1]+1] = meta[j]["_FillValue"]
                rotated[rotated < meta[j]["valid_range"][0]-1] = meta[j]["_FillValue"]
                
                #kvp 
                #values[index][:, :] = np.flipud(rotated) # put into variable
                values[index][0, :, :] = np.flipud(rotated) # put into variable
                #sensors[p_vars].append(np.flipud(rotated)) #add it to complete dataset for LEOGEO calculations

                #time
                end_timer = time.time()
                print("Sensor: ",s_name ," time: ", end_timer - start_timer)
        
        # reset
        indata_temp = []
        inlat_temp = []
        inlon_temp = []
        
    # end individual sensor gridding and adding to nc
    
    # calculate LEOGEO statistics
    # LEOGEO: Mean, STD, NumberOfSensors, SensorWeighting, TotalPixels
    # index of value corresponds to individual sensor
    # leogeo_stats = {"Mean": [[g16], [g17], [viirs], [modis]], "STD": [], "TotalPixels":[]}
    # sensor_idx = {} #key = sensor name, value = sensoridx
    leogeo_calculated_statistics = []
    number_of_sensors_variable = []
    sensor_idx_variable = []
    
    #print("\nLEOGEO STATS\n")
    
    # FilteredQA_550
    # Mean, STD, and TotalPixels stat calculations
    j = 0
    leogeo_index = 0 #index for leogeo statistics array where nc variables are appended to
    leogeo_meta_index = 0
    for p_var in phy_list:
        if p_var in leogeo_stats_arr: #go through each geophysical variable
        
            leogeo_stats = leogeo_stats_arr[p_var] #current dict 
            leogeo_calculated_statistics.append([])
            
            i=0
            # mean, std, pixel statistics
            for statistic in leogeo_stats:
                name = nc_var_name(p_var, "LEOGEO", statistic)
                
                stat_values = []
                if not("Pixels" in statistic) or not("TotalPixels" in statistic):
                   stat_values = calculate_statistic(statistic, leogeo_stats["Mean"])
                else:
                   stat_values = calculate_statistic(statistic, leogeo_stats["TotalPixels"])
               
                #kvp
                #leogeo_calculated_statistics[leogeo_index].append(ds.createVariable(name, np.short,
                #    ('lat', 'lon', ), fill_value = -9999 ,compression='zlib')) #1/29/2023 - added fill value
                leogeo_calculated_statistics[leogeo_index].append(ds.createVariable(name, np.short,
                    ('time', 'lat', 'lon', ), fill_value = -9999 ,compression='zlib')) #1/29/2023 - added fill value

                #print("MAX STAT VALUE for ", statistic, ":", stat_values.max())
                leogeo_long = meta[leogeo_meta_index]["long_name"]
                lfilter = True if ("filtered" in name.lower()) else False
                lname = nc_long_name(p_var, "LEOGEO",statistic,filtered=lfilter )
                affix='AOD_AllQA_550_*_*_Mean'
                if lfilter: affix='AOD_FilteredQA_550_*_*_Mean'

                if statistic == "Mean":
                    leogeo_calculated_statistics[leogeo_index][i].long_name = lname  + \
                    " for the grid. It is the average of individual sensor gridded AODs, i.e., the average of all "+affix+\
                    "  where *_* stands for individual sensor naming."
                else:
                    leogeo_calculated_statistics[leogeo_index][i].long_name = lname  + " for the grid"
                
                # metadata
                leogeo_calculated_statistics[leogeo_index][i].units = "None"  #"degree"#meta[leogeo_meta_index]["units"]
                #print("LEOGEO STAT: ", meta[leogeo_meta_index])
                if not("Pixels" in statistic) or not("TotalPixels" in statistic):
                    leogeo_calculated_statistics[leogeo_index][i].scale_factor = round(meta[leogeo_meta_index]["scale_factor"], 3)
                    leogeo_calculated_statistics[leogeo_index][i].add_offset = meta[leogeo_meta_index]["add_offset"]
                    leogeo_calculated_statistics[leogeo_index][i].valid_range = meta[leogeo_meta_index]["valid_range"]
                else:
                    #Total pixels does nott have scale_factor, add_offset
                    leogeo_calculated_statistics[leogeo_index][i].valid_range = [0, 10000]
                #kvp
                #leogeo_calculated_statistics[leogeo_index][i][:, :] = stat_values
                leogeo_calculated_statistics[leogeo_index][i][0, :, :] = stat_values
                i+=1
            
            
            sensor_order = ["MODIS_Terra", "MODIS_Aqua", "VIIRS_SNPP","VIIRS_NOAA20", "ABI_G16", "ABI_G17", "AHI_H08","AHI_H09"]
            name = nc_var_name(p_var,"LEOGEO", "SensorWeighting")
            sensor_idx_variable.append(ds.createVariable(name, np.short,
                            ('sensor', 'lat', 'lon', ), fill_value=-9999.,compression='zlib'))
            sensor_idx_variable[leogeo_index].long_name = nc_long_name(p_var, "LEOGEO",
                            "SensorWeighting", meta[j]["long_name"]) + " for the grid"
            sensor_idx_variable[leogeo_index].units = "None"
            #kvp
            sensor_idx_variable[leogeo_index].valid_range = [0, 1]
            number_of_sensors = []
            
            i = 0
            for satellite in sensor_order:
                # check to make sure sensor is in there
                if satellite in sensor_idx_arr[p_var]:
                    
                    #add number of sensors together
                    if len(number_of_sensors)==0:
                        number_of_sensors = sensor_idx_arr[p_var][satellite]
                    else:
                        number_of_sensors = number_of_sensors + sensor_idx_arr[p_var][satellite]
                    
                    #save for sensor layer
                    sensor_idx_variable[leogeo_index][i, :, :] = sensor_idx_arr[p_var][satellite]
                # no files for this satellite
                else:
                    pass
                i+=1
            
            name = nc_var_name(p_var, "LEOGEO", "NumberOfSensors")
            #kvp
            #number_of_sensors_variable.append(ds.createVariable(name, np.short, ('lat', 'lon', ),compression='zlib'))
            #number_of_sensors_variable[leogeo_index][:, :] = number_of_sensors
            number_of_sensors_variable.append(ds.createVariable(name, np.short, ('time', 'lat', 'lon', ),compression='zlib'))
            number_of_sensors_variable[leogeo_index][0, :, :] = number_of_sensors

            number_of_sensors_variable[leogeo_index].long_name = nc_long_name(p_var, "LEOGEO",
                                     "NumberOfSensors", meta[j]["long_name"])  + " for the grid"
            number_of_sensors_variable[leogeo_index].units = "None"
            #kvp
            number_of_sensors_variable[leogeo_index].valid_range = [0, len(sensor_order)] # len(full_satellite_list)
            
            leogeo_index+=1    
        j += 1
        leogeo_meta_index += 1
    
    
    return ds
    # output = 1 netcdf with six variables (one for each sensor) 
    # merge = run two loops (lat lon) (xdim, ydim from grid) 
    # merge aod equal to average of six sensors making sure no nan
    # identify which sensors are availabe - sensor ID variable 0 or 1 (unavailable vs available)
    # byte type array, keep n-dimensional for now
    # minimum number of things to merge: in the future

    # future capabilities:
    # minimum number, flagging system, use higher quality data (prefiltered)
