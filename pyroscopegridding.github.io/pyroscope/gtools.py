__author__ = "Sally Zhao"
__copyright__ = "Copyright 2023, Pyroscope"
__credits__ = ["Neil Gutkin", "Jennifer Wei", "Pawan Gupta", "Robert Levy", "Xiaohua Pan", "Zhaohui Zhang"]
__version__ = "1.0.0"
__maintainer__ = "Sally Zhao"
__email__ = "zhaosally0@gmail.com"
__status__ = "Production"
# Gtools
#
# the control for the file. It parses command line arguments and calls to grid_ncf
#

#import other module files
from concurrent.futures import process
import pyroscope.sat_data_input as sat_data_input
import pyroscope.filter_data as filter_data
import pyroscope.gridding as gridding
import pyroscope.grid_ncf as grid_ncf
import pyroscope.fileparser as fileparser # comment for specific file functions
import pyroscope.time_conv as time_conv

import os
import argparse
import numpy as np
import datetime
import time
"""
# python3 gtools.py -r -fn "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/AERDT_L2_VIIRS_SNPP.A2021008.2018.001.2021009072852.nc" -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean"
def read(filename, geo_list, phy_list, verbose = False):

    geo_list = geo_list.split(" ")
    print("Geolocation variables: ", geo_list)
    phy_list = phy_list.split(" ")
    print("Geophysical variables: ", phy_list)
    # open file
    L2FID,GeoID,PhyID =sat_data_input.open_file(filename)

    geo = sat_data_input.read_geo_data(GeoID, var_list=geo_list)
    phy = sat_data_input.read_phy_data(PhyID, GeoID, var_list=phy_list, var_nc=phy_list, var_hdf=phy_list)
    lat = geo['data'][0][:,:]
    lon = geo['data'][1][:,:]
    aod = phy['data'][0][:,:]

    row,col=np.array(aod).shape
    

     # valid data points
    lat1=lat[aod>0]
    lon1=lon[aod>0]
    aod1=aod[aod>0]

    L2FID.close()

    if verbose:
        print("Gridding: " + str(filename) + "\nGeolocation:" + str(geo_list)+ "\nGeophysical:" + str(phy_list))

        print("dims:", row, col)
        #PRINTS
        for i in range(10): #(aod1.size):
            print('(lat,lon)=({0:.2f}, {1:.2f}) aod= {2:.3f}'.\
                format(lat1[i], lon1[i], aod1[i]))
    
    return lat1, lon1, aod1

# python3 gtools.py -f -fn "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/AERDT_L2_VIIRS_SNPP.A2021008.2018.001.2021009072852.nc" -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean"
def filter(filename, geo_list, phy_list, verbose = False):
    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)

    # open file
    L2FID,GeoID,PhyID = sat_data_input.open_file(filename)
    #filter
    lat,lon,phy_vars, metadata = filter_data.filter_data(GeoID, PhyID, geo_list, phy_list, phy_nc=phy_list, phy_hdf=phy_list)
    # close file
    L2FID.close()
    
    if verbose:
        print("Latitude: ", lat)
        print("Latitude: ", lon)
        print("Geophysical: ", phy_vars)

    return lat,lon,phy_vars, metadata

def grid(filename, geo_list, phy_list, verbose = False):
    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)

    L2FID,GeoID,PhyID = sat_data_input.open_file(filename)

    lat,lon,phy_vars, metadata = filter_data.filter_data(GeoID, PhyID, geo_list, phy_list, phy_nc=phy_list, phy_hdf=phy_list)

    L2FID.close()

    #grid parameters
    limit = [min(lat[0]), max(lat[0]), min(lon[0]),  max(lon[0])]
    gsize = 0.25
    indata = phy_vars[0]
    inlat=lat[0]
    inlon = lon[0]

    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,gsize,indata,inlat,inlon)

    if verbose:
        print(avgtau)

    return avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau
"""
"""
def netCDF_single(filename, geo_list, phy_list, output):
    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)

    #L2FID,GeoID,PhyID = sat_data_input.open_file(filename)

    #,lon,phy_vars, metadata = filter_data.filter_data(GeoID, PhyID, geo_list, phy_list, phy_nc=phy_list, phy_hdf=phy_list)

    #L2FID.close()

    #grid parameters
    limit = [min(lat[0]), max(lat[0]), min(lon[0]),  max(lon[0])]
    gsize = 0.25
    indata = phy_vars[0]
    inlat=lat[0]
    inlon = lon[0]

    print("Results at: " + output)

    grid_ncf.gridNC_single(limit, gsize, indata, inlat, inlon, output)
"""
#gridNC_time(limit, gsize, indata, inlat, inlon, 1, phy_list, geo_list, sat_files[:10], fpath + "viirs_time_compv1.nc")
def netCDF_multi(filelist, geo_list, phy_list, output):
    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)
    filelist = fileparser.read_file_sat_data(filelist)

    print(filelist)

    lat,lon,phy_vars = gridding.multi_sensor_grid_data(filelist, phy_list, geo_list)

    print(lat)
    print(lon)
    print(phy_vars)

    #grid parameters
    limit = [min(lat), max(lat), min(lon),  max(lon)]
    gsize = 0.25
    indata = phy_vars
    inlat=lat
    inlon = lon

    print("Results at: " + output)

    #grid_ncf.gridNC_time(limit, gsize, indata, inlat, inlon, 1, phy_list, geo_list, filelist, output)
    
    #grid_ncf.gridNC_time(limit, gsize, indata, inlat, inlon, output)


    
def netCDF_multi_time(filelist, gsize, geo_list, phy_list, output, time_interval, time_start, time_end):
    start_timer = time.time()

    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)
    filelist = fileparser.read_file_sat_data(filelist)

    start = time_conv.to_datetime(time_start)
    end = time_conv.to_datetime(time_end)

    split_files = time_conv.split_filetimes(filelist, start, end, int(time_interval))

    lat,lon,phy_vars = gridding.multi_sensor_grid_data(filelist, phy_list, geo_list)
    #grid parameters
    limit = [min(lat), max(lat), min(lon),  max(lon)]
    gsize = float(gsize)
    indata = phy_vars
    inlat=lat
    inlon = lon

    #grid_ncf.grid_nc_mult_files_time(limit, gsize, indata, inlat, inlon, phy_list, geo_list, split_files, output)

    end_timer = time.time()
    print("Full execution time: ", end_timer - start_timer)

# output - folder location
def netCDF_sensor_statistics(filelist, limit, gsize, geo_list, phy_list, outputfolder, outputname, time_interval, time_start, time_end):
    start_timer = time.time()
    
    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)
    limit = [float(item) for item in limit.split(" ")]
    filelist = fileparser.read_file_sat_data(filelist)

    start = time_conv.to_datetime(time_start)
    end = time_conv.to_datetime(time_end)

    split_files = time_conv.split_filetimes(filelist, start, end, int(time_interval))

    gsize = float(gsize)

    curr = start
    for f in split_files:
        name = outputfolder + outputname + str(curr).replace(":", "-")+".nc"

        # grid files
        #grid_ncf.grid_nc_single_statistics(limit, gsize, geo_list, phy_list, f, name)
        curr = curr + datetime.timedelta(minutes = int(time_interval))

    end_timer = time.time()
    print("Full execution time: ", end_timer - start_timer)

# output - folder location
def netCDF_sensor_statistics_split(filelist, limit, gsize, geo_list, phy_list, outputfolder, outputname, time_interval, time_start, time_end):

    start_timer = time.time()

    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)
    limit = [float(item) for item in limit.split(" ")]
    filelist = fileparser.read_file_sat_data(filelist)

    start = time_conv.to_datetime(time_start)
    end = time_conv.to_datetime(time_end)

    # [[filenames for a time interval] , [], []]
    split_files = time_conv.split_filetimes(filelist, start, end, int(time_interval))
    
    # [{sensor: array of filenames}, {}, {}, {}, {}]
    # code should record bugs - always six sensors available
    # code should report that modis is not available
    split_files = time_conv.split_filenames(split_files)

    gsize = float(gsize)

    curr = start
    for f in split_files:
        name = outputfolder + outputname + str(curr).replace(":", "-")+".nc"

        # grid files
        #grid_ncf.grid_nc_statistics(limit, gsize, geo_list, phy_list, f, name)
        curr = curr + datetime.timedelta(minutes = int(time_interval))

    end_timer = time.time()
    print("Full execution time: ", end_timer - start_timer)

# output - folder location
def netCDF_sensor_statistics_id(filelist, limit, gsize, geo_list, phy_list, outputfolder, outputname, time_interval, time_start, time_end):

    start_timer = time.time()

    geo_list = geo_list.split(" ")
    print(geo_list)
    phy_list = phy_list.split(" ")
    print(phy_list)
    limit = [float(item) for item in limit.split(" ")]
    filelist = fileparser.read_file_sat_data(filelist)

    start = time_conv.to_datetime(time_start)
    end = time_conv.to_datetime(time_end)

    # [[filenames for a time interval] , [], []]
    split_files = time_conv.split_filetimes(filelist, start, end, int(time_interval))
    
    # [{sensor: array of filenames}, {}, {}, {}, {}]
    # code should record bugs - always six sensors available
    # code should report that modis is not available
    split_files = time_conv.split_filenames(split_files)

    gsize = float(gsize)

    curr = start
    for f in split_files:
        name = outputfolder + outputname + str(curr).replace(":", "-")+".nc"
        #print(f)
        # print(split_files)
        
        # grid files
        #grid_nc_sensor_statistics_metadata
        #grid_ncf.grid_nc_sensor_statistics_metadata(limit, gsize, geo_list, phy_list, f, name, phy_nc=phy_list, phy_hdf=phy_list)
        curr = curr + datetime.timedelta(minutes = int(time_interval))

    end_timer = time.time()
    print("Full execution time: ", end_timer - start_timer)

# output - folder location
def netCDF_yaml_config(file_name):

    grid_settings, variables, file_io = fileparser.read_config(file_name)

    start_timer = time.time()

    geo_list = variables["geo_var"]
    phy_list = variables["phy_var"]
    phy_nc = variables["phy_var_nc"]
    phy_hdf = variables["phy_var_hdf"]
    limit = grid_settings["limit"]
    pixel_range = variables["pixel_range"]
    
    # list of file paths for specific files to read
    if file_io["file_directory_folder"] == "NA":
        if file_io["file_location_folder"] == "NA":
            filelist = fileparser.read_file_sat_data(file_io["file_location_file"])
        else:
            filelist = fileparser.read_folder_sat_data(file_io["file_location_folder"])
    else:
        filelist = fileparser.read_directory_sat_data(file_io["file_directory_folder"])
        
        
    #static file for Land_Sea_Mask and Topographic_Altitude
    static_file = file_io["static_file"]
    print("static file:", static_file)

    time_start = grid_settings["time_start"]
    time_end = grid_settings["time_end"]
    start = time_conv.to_datetime(time_start)
    end = time_conv.to_datetime(time_end)

    # [[filenames for a time interval] , [], []]
    time_interval = grid_settings["time_interval"]
    split_files = time_conv.split_filetimes(filelist, start, end, int(time_interval))
    
    # [{sensor: array of filenames}, {}, {}, {}, {}]
    # code should record bugs - always six sensors available
    # code should report that modis is not available
    split_files = time_conv.split_filenames(split_files)

    gsize = float(grid_settings["gridsize"])

    curr = start
    outputfolder = file_io["output_location"]
    outputname = file_io["output_name"]
    
    count= 0
    for f in split_files:
        #name = outputfolder + outputname + str(curr).replace(":", "-")+".nc"
        #XAERDT_L3_MEASURES_QD_HH.YYYYMMDD.HHMM.V0.ProcessingDate.nc
        time_string = str(curr).replace(":", "").replace("-", "").replace(" ", "")
        time_string = time_string[:8] + '.' + time_string[8:-2]
        processing_date = time_conv.sys_time()
        name = outputfolder + "XAERDT_L3_MEASURES_QD_HH." + time_string+ ".V0." +processing_date+ ".nc"
        
        #print(f) f = dict (key: sensor, values: array of filelocs)
        # print(split_files)
        
        #print("\n\nNAME:", name)
        #print("\n\TIME STRING:", curr)
        #print("\n\PROCESS DATE:", time_interval)
        
        # grid files
        #grid_nc_sensor_statistics_metadata
        
        ds = grid_ncf.grid_nc_sensor_statistics_metadata(limit, gsize, geo_list, phy_list, f, name, curr, time_interval, phy_nc, phy_hdf, static_file, pixel_range)
        ds.close()
        curr = curr + datetime.timedelta(minutes = int(time_interval))

    end_timer = time.time()
    print("Full execution time: ", end_timer - start_timer)
    
# output - folder location
def netCDF_yaml_config_time(file_name, time_start, time_end):

    grid_settings, variables, file_io = fileparser.read_config(file_name)

    start_timer = time.time()

    geo_list = variables["geo_var"]
    phy_list = variables["phy_var"]
    phy_nc = variables["phy_var_nc"]
    phy_hdf = variables["phy_var_hdf"]
    limit = grid_settings["limit"]
    pixel_range = variables["pixel_range"]
    
    # list of file paths for specific files to read
    if file_io["file_directory_folder"] == "NA":
        if file_io["file_location_folder"] == "NA":
            filelist = fileparser.read_file_sat_data(file_io["file_location_file"])
        else:
            filelist = fileparser.read_folder_sat_data(file_io["file_location_folder"])
    else:
        filelist = fileparser.read_directory_sat_data(file_io["file_directory_folder"])
        
    print("FILELIST: ", filelist)
    #static file for Land_Sea_Mask and Topographic_Altitude
    static_file = file_io["static_file"]
    print("static file:", static_file)

    #time_start = grid_settings["time_start"]
    #time_end = grid_settings["time_end"]
    start = time_conv.to_datetime(time_start)
    end = time_conv.to_datetime(time_end)
    

    # [[filenames for a time interval] , [], []]
    time_interval = grid_settings["time_interval"]
    split_files = time_conv.split_filetimes(filelist, start, end, int(time_interval))
    
    # [{sensor: array of filenames}, {}, {}, {}, {}]
    # code should record bugs - always six sensors available
    # code should report that modis is not available
    split_files = time_conv.split_filenames(split_files)

    gsize = float(grid_settings["gridsize"])

    curr = start
    outputfolder = file_io["output_location"]
    outputname = file_io["output_name"]
    
    count= 0
    for f in split_files:
        #print("FILELIST:",f)
        #f = {'AERDT_L2_ABI_G16': ['/mnt/c/Users/bobgr/Desktop/NASA Spring 2023/Gridtools Package (Code, README, inputs, outputs, examples, verification)/LAADSemulation/g16/XAERDT_L2_ABI_G16.A2019208.0000.001.2022260004939.nc', '/mnt/c/Users/bobgr/Desktop/NASA Spring 2023/Gridtools Package (Code, README, inputs, outputs, examples, verification)/LAADSemulation/g16/XAERDT_L2_ABI_G16.A2019208.0010.001.2022260004939.nc.nc', '/mnt/c/Users/bobgr/Desktop/NASA Spring 2023/Gridtools Package (Code, README, inputs, outputs, examples, verification)/LAADSemulation/g16/XAERDT_L2_ABI_G16.A2019208.0020.001.2022260004940.nc']}
        #name = outputfolder + outputname + str(curr).replace(":", "-")+".nc"
        #XAERDT_L3_MEASURES_QD_HH.YYYYMMDD.HHMM.V0.ProcessingDate.nc
        time_string = str(curr).replace(":", "").replace("-", "").replace(" ", "")
        time_string = time_string[:8] + '.' + time_string[8:-2]
        processing_date = time_conv.sys_time()
        name = outputfolder + "XAERDT_L3_MEASURES_QD_HH." + time_string+ ".V0." +processing_date+ ".nc"
        
        ds = grid_ncf.grid_nc_sensor_statistics_metadata(limit, gsize, geo_list, phy_list, f, name, curr, time_interval, phy_nc, phy_hdf, static_file, pixel_range)
        ds.close()
        curr = curr + datetime.timedelta(minutes = int(time_interval))

    end_timer = time.time()
    print("Full execution time: ", end_timer - start_timer)

    

def execute_file(filename):
    filelist, gsize, time_interval, start, end, output, geo_list, phy_list = fileparser.read_in_commands(filename)
    lat,lon,phy_vars = gridding.multi_sensor_grid_data(filelist, phy_list, geo_list)
    #grid parameters
    limit = [min(lat), max(lat), min(lon),  max(lon)]
    gsize = float(gsize)
    indata = phy_vars
    inlat=lat
    inlon = lon

    split_files = time_conv.split_filetimes(filelist, start, end, int(time_interval))

    #grid_ncf.grid_nc_mult_files_time(limit, gsize, indata, inlat, inlon, phy_list, geo_list, split_files, output)



def control():
    # parser object
    parser = argparse.ArgumentParser()

    # reads in commands as mutually exclusive group
    cmd_group = parser.add_mutually_exclusive_group()
    cmd_group.add_argument("-r", "--read", action="store_true", help="Reads raw data given file")
    cmd_group.add_argument("-f", "--filter", action="store_true", help="Reads and filters data given file")
    cmd_group.add_argument("-g", "--grid", action="store_true", help="Reads, filters, and grids data given file")
    cmd_group.add_argument("-ns", "--netCDF_single", action="store_true", help="Saves gridded as netCDF given single file")
    cmd_group.add_argument("-nm", "--netCDF_multi", action="store_true", help="Saves gridded as netCDF given multiple files. Fuses data regardless of time.")
    cmd_group.add_argument("-nmt", "--netCDF_multi_time", action="store_true", help="Saves gridded as netCDF given multiple files. \nFuses data according to time interval.")
    cmd_group.add_argument("-exec", "--execute_file", action="store_true", help="Reads commands from textfile")
    cmd_group.add_argument("-ss", "--sensor_statistics", action="store_true", help="Reads sensors and reports statistics and individual gridded data.")
    cmd_group.add_argument("-sss", "--sensor_statistics_split", action="store_true", help="Reads sensors and reports statistics and gridded data based on satellite categorization.")
    cmd_group.add_argument("-ssi", "--sensor_statistics_split_id", action="store_true", help="Reads sensors and reports statistics and gridded data based on satellite categorization.")
    cmd_group.add_argument("-cfg", "--config", action="store_true", help="Grids according to yaml config file")
    cmd_group.add_argument("-cfgtime", "--configtime", action="store_true", help="Grids according to yaml config file")


    # verbose, quiet
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument("-v", "--verbose", action="store_true", help="Outputs verbose results")
    output_group.add_argument("-q", "--quiet", action="store_true", help="Outputs quiet results")

    # parameters
    parser.add_argument("-fn", "--filename", help="Filename reading in")
    parser.add_argument("-fl", "--filelist", help="Location of text file containing filepaths reading in")
    parser.add_argument("-gl", "--geo_list", help="Geolocation variables (if modis, [Latitiude, Longitude])")
    parser.add_argument("-gp", "--geo_phys", help="Geophysical variables")
    parser.add_argument("-gs", "--gsize", help="Gridding size")
    parser.add_argument("-o", "--output", help = "Output result to file")
    parser.add_argument("-on", "--output_name", help = "Output naming convention for multiple files")
    parser.add_argument("-l", "--limit", help = "Boundary box for gridding")

    # time variables
    parser.add_argument("-ti", "--time_interval", help="Time Interval (minutes)")
    parser.add_argument("-ts", "--time_start", help="Time start ('yyyy/mm/dd/hr/mm')")
    parser.add_argument("-te", "--time_end", help="Time end ('yyyy/mm/dd/hr/mm')")

    

    args = parser.parse_args()

    # decides what function to run
    if args.read:
        pass
        #python3 gtools.py -r -fn "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleInputs 0000-0059 01-01-2020/AERDT_L2_ABI_G16.A2020001.0000.001.2022003073056.nc" -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean"
        #read(str(args.filename), str(args.geo_list),str(args.geo_phys), args.verbose)
    elif args.filter:
        pass
        #python3 gtools.py -f -fn "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleInputs 0000-0059 01-01-2020/AERDT_L2_ABI_G16.A2020001.0000.001.2022003073056.nc" -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean"
        #filter(str(args.filename), str(args.geo_list),str(args.geo_phys), args.verbose)
    elif args.grid:
        pass
        #python3 gtools.py -g -fn "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleInputs 0000-0059 01-01-2020/AERDT_L2_ABI_G16.A2020001.0000.001.2022003073056.nc" -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean"
        #grid(str(args.filename), str(args.geo_list),str(args.geo_phys), args.verbose)
    elif args.netCDF_single:
        pass
        #python3 gtools.py -ns -fn "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleInputs 0000-0059 01-01-2020/AERDT_L2_ABI_G16.A2020001.0000.001.2022003073056.nc" -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean" -o "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleOutputs 0000-0059 01-01-2020/single_ns_result.nc"
        #netCDF_single(str(args.filename), str(args.geo_list),str(args.geo_phys), args.output)
    elif args.netCDF_multi: # fuses these together without regard to time
        #python3 gtools.py -nm -fl "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/satellite_paths.txt"  -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean" -o "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleOutputs 0000-0059 01-01-2020/multi_ns_result.nc"
        netCDF_multi(str(args.filelist), str(args.geo_list),str(args.geo_phys), args.output)
    elif args.netCDF_multi_time: # fuses these together with regard to time
        #python3 gtools.py -nmt -fl "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/satellite_paths.txt" -gs 0.25 -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean" -o "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleOutputs 0000-0059 01-01-2020/multi_ns_time_result.nc" -ti 30 -ts "2020/01/01/00/00" -te "2020/01/01/02/00"
        netCDF_multi_time(str(args.filelist),  str(args.gsize), str(args.geo_list),str(args.geo_phys), args.output, 
            str(args.time_interval), str(args.time_start), str(args.time_end))
    elif args.sensor_statistics:
        #python3 gtools.py -ss -fl "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/satellite_paths.txt" -l "-90 90 -180 180" -gs 0.25 -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean" -o "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleOutputs 0000-0059 01-01-2020/" -on "STATS" -ti 30 -ts "2020/01/01/00/00" -te "2020/01/01/01/00" 
        netCDF_sensor_statistics(str(args.filelist), str(args.limit), str(args.gsize), 
                str(args.geo_list), str(args.geo_phys), args.output, args.output_name, 
                str(args.time_interval), str(args.time_start), str(args.time_end))
    elif args.sensor_statistics_split: 
        # python3 gtools.py -sss -fl "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/satellite_paths.txt" -l "-90 90 -180 180" -gs 0.25 -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean" -o "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleOutputs 0000-0059 01-01-2020/" -on "STATS_SPLIT" -ti 30 -ts "2020/01/01/00/00" -te "2020/01/01/01/00"
        netCDF_sensor_statistics_split(str(args.filelist), str(args.limit), str(args.gsize), 
                str(args.geo_list), str(args.geo_phys), args.output, args.output_name, 
                str(args.time_interval), str(args.time_start), str(args.time_end))
    elif args.sensor_statistics_split_id: 
        #python3 gtools.py -ssi -fl "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/satellite_paths.txt" -l "-90 90 -180 180" -gs 0.25 -gl "latitude longitude" -gp "Optical_Depth_Land_And_Ocean" -o "/mnt/c/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleOutputs 0000-0059 01-01-2020/" -on "STATS_SENSOR_SPLIT" -ti 30 -ts "2020/01/01/00/00" -te "2020/01/01/01/00"
        netCDF_sensor_statistics_id(str(args.filelist), str(args.limit), str(args.gsize), 
                str(args.geo_list), str(args.geo_phys), args.output, args.output_name, 
                str(args.time_interval), str(args.time_start), str(args.time_end))
    elif args.config:
        # yaml config 
        # python gtools.py -cfg -fn "C:\Users\swzhao\Desktop\repos\gridtools\config.yml"
        # "C:\Users\bobgr\Desktop\Spring 2022 NASA\Gridtools Package (Code, README, inputs, outputs, examples, verification)\gridtools\config.yml"
        # "C:/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/gridtools/config.yml"
        netCDF_yaml_config(str(args.filename))
    elif args.configtime:
        # yaml config 
        # python gtools.py -cfg -fn "C:\Users\swzhao\Desktop\repos\gridtools\config.yml"
        # "C:\Users\bobgr\Desktop\Spring 2022 NASA\Gridtools Package (Code, README, inputs, outputs, examples, verification)\gridtools\config.yml"
        # "C:/Users/bobgr/Desktop/Spring 2022 NASA/Gridtools Package (Code, README, inputs, outputs, examples, verification)/gridtools/config.yml"
        netCDF_yaml_config_time(str(args.filename), str(args.time_start), str(args.time_end))
    else:
        print("No command given")
        

if __name__ == '__main__':
    control()
    


