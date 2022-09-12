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

import netCDF4
import numpy as np
import numpy.ma as ma
import time

import gridding
import sat_data_input
import filter_data

#filename - 'gridNet.nc'
# takes representation of a single file
# to make geo2D variable - add _CoordinateAxisType = "Lat";
def gridNC_single(limit, gsize, indata, inlat, inlon, filename):

    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    #lat_dim= len(np.arange(minlat, maxlat, gsize)) #int(1+round(((maxlat-minlat)/gsize)))
    #lon_dim= len(np.arange(minlon, maxlon, gsize)) # int(1+round(((maxlon-minlon)/gsize)))
    lat_dim = int(1+round(((maxlat-minlat)/gsize)))
    lon_dim = int(1+round(((maxlon-minlon)/gsize)))
    
    #set target netcdf file
    ds = netCDF4.Dataset(filename, 'w',format='NETCDF4')
    
    # Create time, lat, lon dimensions
    time = ds.createDimension('time', None)
    lat = ds.createDimension('lat', lat_dim) #latitude is fatitude (y)
    lon = ds.createDimension('lon', lon_dim) # x
    
    #set time variable to filename parse
    times = ds.createVariable('time', 'f4', ('time',))
    
    #set lats and lons  (or rather y and x)
    lats = ds.createVariable('lat', 'f4', ('lat',)) #latitude is fatitude
    lons = ds.createVariable('lon', 'f4', ('lon',))
    
    value = ds.createVariable('value', 'f4', ('time', 'lat', 'lon', ))
    value.units = 'unknown'
    
    #bounds for lats and lons
    lats[:] = np.arange(0,lat_dim)*gsize+minlat
    lons[:] = np.arange(0,lon_dim)*gsize+minlon
    
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,gsize,indata,inlat,inlon)
    
    rotated = list(zip(*mintau))[::-1]
    value[0, :, :] = (np.flipud(rotated))
    
    #close file
    ds.close()
    
# filename - 'gridNet.nc' (this is what it saves to)
# sat_files - list of satellite files to pull data from
def gridNC_time(limit, gsize, indata, inlat, inlon, time_interval, phy_list, geo_list, sat_files, filename):

    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    #lat_dim= len(np.arange(minlat, maxlat, gsize)) #int(1+round(((maxlat-minlat)/gsize)))
    #lon_dim= len(np.arange(minlon, maxlon, gsize)) # int(1+round(((maxlon-minlon)/gsize)))
    lat_dim = int(1+round(((maxlat-minlat)/gsize)))
    lon_dim = int(1+round(((maxlon-minlon)/gsize)))
    
    #set target netcdf file
    ds = netCDF4.Dataset(filename, 'w',format='NETCDF4')
    
    # Create time, lat, lon dimensions
    time = ds.createDimension('time', None)
    lat = ds.createDimension('lat', lat_dim) #latitude is fatitude (y)
    lon = ds.createDimension('lon', lon_dim) # x
    
    #set time variable to filename parse
    times = ds.createVariable('time', 'f4', ('time',))
    
    #set lats and lons  (or rather y and x)
    lats = ds.createVariable('lat', 'f4', ('lat',)) #latitude is fatitude
    lons = ds.createVariable('lon', 'f4', ('lon',))
    
    value = ds.createVariable('value', 'f4', ('time', 'lat', 'lon', ))
    value.units = 'unknown'
    
    #bounds for lats and lons
    lats[:] = np.arange(0,lat_dim)*gsize+minlat
    lons[:] = np.arange(0,lon_dim)*gsize+minlon
    
    #run through and assign gridded data to corresponding time interval    
    for i in range(0,len(sat_files),time_interval):
        satellites = sat_files[i: i+time_interval]
        print(satellites)
        avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau =  gridding.multi_sensor_grid(satellites, gsize, limit, phy_list, geo_list)
    
        rotated = list(zip(*mintau))[::-1]
        value[i, :, :] = (np.flipud(rotated))
        
    
    #close file
    ds.close()

# filename - 'gridNet.nc' (this is what it saves to)
# sat_files - list of satellite files to pull data from
def grid_nc_mult_files_time(limit, gsize, indata, inlat, inlon, phy_list, geo_list, filelist, filename):

    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    lat_dim = int(1+round(((maxlat-minlat)/gsize)))
    lon_dim = int(1+round(((maxlon-minlon)/gsize)))
    
    #set target netcdf file
    ds = netCDF4.Dataset(filename, 'w',format='NETCDF4')
    
    # Create time, lat, lon dimensions
    time = ds.createDimension('time', None)
    lat = ds.createDimension('lat', lat_dim) #latitude is fatitude (y)
    lon = ds.createDimension('lon', lon_dim) # x
    
    #set time variable to filename parse
    times = ds.createVariable('time', 'f4', ('time',))
    
    #set lats and lons  (or rather y and x)
    lats = ds.createVariable('lat', 'f4', ('lat',)) #latitude is fatitude
    lons = ds.createVariable('lon', 'f4', ('lon',))
    
    value = ds.createVariable('value', 'f4', ('time', 'lat', 'lon', ))
    value.units = 'unknown'
    
    #bounds for lats and lons
    lats[:] = np.arange(0,lat_dim)*gsize+minlat
    lons[:] = np.arange(0,lon_dim)*gsize+minlon
    
    #run through and assign gridded data to corresponding time interval    
    for i in range(0,len(filelist)):
        avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau =  gridding.multi_sensor_grid(filelist[i], gsize, limit, phy_list, geo_list)
    
        rotated = list(zip(*mintau))[::-1]
        value[i, :, :] = (np.flipud(rotated))
        
    
    #close file
    ds.close()

# filename - 'gridNet.nc' (this is what it saves to)
# sat_files - list of satellite files to pull data from
# saves all satellite sensor gridded data as separate variables
#def grid_nc_single_statistics(limit, gsize, inlat, inlon, geo_list, phy_list, filelist, filename):
def grid_nc_single_statistics(limit, gsize, geo_list, phy_list, filelist, filename):

    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    lat_dim = int(1+round(((maxlat-minlat)/gsize)))
    lon_dim = int(1+round(((maxlon-minlon)/gsize)))
    
    #set target netcdf file
    ds = netCDF4.Dataset(filename, 'w',format='NETCDF4')
    
    # Create time, lat, lon dimensions
    timed = ds.createDimension('time', None)
    lat = ds.createDimension('lat', lat_dim) #latitude is fatitude (y)
    lon = ds.createDimension('lon', lon_dim) # x
    
    #set time variable to filename parse
    times = ds.createVariable('time', 'f4', ('time',))
    
    #set lats and lons  (or rather y and x)
    lats = ds.createVariable('lat', 'f4', ('lat',)) #latitude is fatitude
    lons = ds.createVariable('lon', 'f4', ('lon',))

    #bounds for lats and lons
    lats[:] = np.arange(0,lat_dim)*gsize+minlat
    lons[:] = np.arange(0,lon_dim)*gsize+minlon
    
    # unknown number of sensors, save values to list
    values = []
    # instantiate empty arrays for summary statistics
    indata_stat = []
    inlat_stat = []
    inlon_stat = []
    
    #run through and assign gridded data to single sensor variable
    for i in range(0,len(filelist)):
        start_timer = time.time()
        name = str(filelist[i].split("/")[-1])
        values.append(ds.createVariable(name, 'f4', ('time', 'lat', 'lon', )))
        values[i].units = 'unknown'

        #print(filelist[i])
        
        L2FID,GeoID,PhyID = sat_data_input.open_file(filelist[i])
        # check for hdf files 
        if (filelist[i].split(".")[-1] == "hdf"):
            geo_list = ['Latitude', 'Longitude']
        else:
            geo_list = ['latitude', 'longitude']

        lat,lon,phy_vars, metadata = filter_data.filter_data(GeoID, PhyID, geo_list, phy_list, phy_nc=phy_list, phy_hdf=phy_list)


        L2FID.close()

        #summary statistics appending
        if len(indata_stat) == 0:
            indata_stat = phy_vars[0]
            inlat_stat = lat[0]
            inlon_stat = lon[0]
        else:
            inlat_stat = ma.concatenate([inlat_stat, lat[0]])
            # dictionary later on
            indata_stat = ma.concatenate([indata_stat, phy_vars[0]]) 
            inlon_stat = ma.concatenate([inlon_stat, lon[0]])

        #grid parameters
        #limit = [min(lat[0]), max(lat[0]), min(lon[0]),  max(lon[0])]
        gsize = 0.25
        indata = phy_vars[0]
        inlat=lat[0]
        inlon = lon[0]

        #dont need this
        avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,gsize,indata,inlat,inlon)

        rotated = list(zip(*avgtau))[::-1]
        #print(np.flipud(rotated))
        #print(i, ": ", filelist[i])
        #values[i][0, :, :] = np.nan_to_num(np.flipud(rotated))
        values[i][0, :, :] = np.flipud(rotated)

        end_timer = time.time()
        print("Sensor: ",name ," time: ", end_timer - start_timer)

    start_timer = time.time()

    # put in gridded statistics
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,gsize,indata_stat,inlat_stat,inlon_stat)

    avg = ds.createVariable("merge_avgtau", 'f4', ('time', 'lat', 'lon', ))
    avg.units = 'unknown' # pull from metadata
    avg[0, :, :] = np.flipud(list(zip(*avgtau))[::-1]) # make sure this is correct

    std = ds.createVariable("merge_stdtau", 'f4', ('time', 'lat', 'lon', ))
    std.units = 'unknown'
    std[0, :, :] = np.flipud(list(zip(*stdtau))[::-1])

    grdlatvar = ds.createVariable("merge_grdlat", 'f4', ('time', 'lat', 'lon', ))
    grdlatvar.units = 'unknown'
    grdlatvar[0, :, :] = np.flipud(list(zip(*grdlat))[::-1])

    grdlonvar = ds.createVariable("merge_grdlon", 'f4', ('time', 'lat', 'lon', ))
    grdlonvar.units = 'unknown'
    grdlonvar[0, :, :] = np.flipud(list(zip(*grdlon))[::-1])

    min = ds.createVariable("merge_mintau", 'f4', ('time', 'lat', 'lon', ))
    min.units = 'unknown'
    min[0, :, :] = np.flipud(list(zip(*mintau))[::-1])

    max = ds.createVariable("merge_maxtau", 'f4', ('time', 'lat', 'lon', ))
    max.units = 'unknown'
    max[0, :, :] = np.flipud(list(zip(*maxtau))[::-1])

    countvar = ds.createVariable("merge_count", 'f4', ('time', 'lat', 'lon', ))
    countvar.units = 'unknown'
    countvar[0, :, :] = np.flipud(list(zip(*count))[::-1])

    sum = ds.createVariable("merge_sumtau", 'f4', ('time', 'lat', 'lon', ))
    sum.units = 'unknown'
    sum[0, :, :] = np.flipud(list(zip(*sumtau))[::-1])

    end_timer = time.time()
    print("Statistics summary time: ", end_timer - start_timer)
    
    #close file
    ds.close()


# filename - 'gridNet.nc' (this is what it saves to)
# sat_files - list of satellite files to pull data from
# saves all satellite sensor gridded data as separate variables
#def grid_nc_single_statistics(limit, gsize, inlat, inlon, geo_list, phy_list, filelist, filename):
def grid_nc_statistics(limit, gsize, geo_list, phy_list, filelist, filename):
    #start_timer = time.time()

    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    lat_dim = int(1+round(((maxlat-minlat)/gsize)))
    lon_dim = int(1+round(((maxlon-minlon)/gsize)))
    
    #set target netcdf file
    ds = netCDF4.Dataset(filename, 'w',format='NETCDF4')
    
    # Create time, lat, lon dimensions
    timed = ds.createDimension('time', None)
    lat = ds.createDimension('lat', lat_dim) #latitude is fatitude (y)
    lon = ds.createDimension('lon', lon_dim) # x
    
    #set time variable to filename parse
    times = ds.createVariable('time', 'f4', ('time',))
    
    #set lats and lons  (or rather y and x)
    lats = ds.createVariable('lat', 'f4', ('lat',)) #latitude is fatitude
    lons = ds.createVariable('lon', 'f4', ('lon',))

    #bounds for lats and lons
    lats[:] = np.arange(0,lat_dim)*gsize+minlat
    lons[:] = np.arange(0,lon_dim)*gsize+minlon
    
    # unknown number of sensors, save values to list
    values = []

    # instantiate empty arrays for summary statistics
    indata_stat = []
    inlat_stat = []
    inlon_stat = []
    indata_temp = []
    inlat_temp = []
    inlon_temp = []
    
    #run through and assign gridded data to sensor variable
    # filenames are in the form:
    # [[time0: [sat1], [sat2], ... ], ...[timek: [sat1], [sat2], ....]] 
    # we assume we receive a signle dict for a time interval though

    for i, s_name in enumerate(filelist.keys()): 
        start_timer = time.time()
        name = str(s_name)
        
        values.append(ds.createVariable(name, 'f4', ('time', 'lat', 'lon', )))
        values[i].units = 'unknown'

        for s in filelist.get(s_name):
            
            L2FID,GeoID,PhyID = sat_data_input.open_file(str(s))

            # check for hdf files 
            if (str(s).split(".")[-1] == "hdf"):
                geo_list = ['Latitude', 'Longitude']
            else:
                geo_list = ['latitude', 'longitude']

            lat,lon,phy_vars, metadata = filter_data.filter_data(GeoID, PhyID, geo_list, phy_list, phy_nc=phy_list, phy_hdf=phy_list)

            L2FID.close()

            # statistics appending for this sensor
            if len(indata_temp) == 0:
                indata_temp = phy_vars[0]
                inlat_temp = lat[0]
                inlon_temp = lon[0]
            else:
                inlat_temp = ma.concatenate([inlat_temp, lat[0]])
                indata_temp = ma.concatenate([indata_temp, phy_vars[0]]) 
                inlon_temp = ma.concatenate([inlon_temp, lon[0]])

            #summary statistics appending
            if len(indata_stat) == 0:
                indata_stat = phy_vars[0]
                inlat_stat = lat[0]
                inlon_stat = lon[0]
            else:
                inlat_stat = ma.concatenate([inlat_stat, lat[0]])
                indata_stat = ma.concatenate([indata_stat, phy_vars[0]]) 
                inlon_stat = ma.concatenate([inlon_stat, lon[0]])

        #grid parameters
        avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,float(gsize),indata_temp,inlat_temp,inlon_temp)

        # reset
        indata_temp = []
        inlat_temp = []
        inlon_temp = []

        rotated = list(zip(*avgtau))[::-1]

        values[i][0, :, :] = np.flipud(rotated)

        end_timer = time.time()
        print("Sensor: ",s_name ," time: ", end_timer - start_timer)

    # put in gridded statistics
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,float(gsize),indata_stat,inlat_stat,inlon_stat)

    avg = ds.createVariable("merge_avgtau", 'f4', ('time', 'lat', 'lon', ))
    avg.units = 'unknown'
    avg[0, :, :] = np.flipud(list(zip(*avgtau))[::-1])

    std = ds.createVariable("merge_stdtau", 'f4', ('time', 'lat', 'lon', ))
    std.units = 'unknown'
    std[0, :, :] = np.flipud(list(zip(*stdtau))[::-1])

    grdlatvar = ds.createVariable("merge_grdlat", 'f4', ('time', 'lat', 'lon', ))
    grdlatvar.units = 'unknown'
    grdlatvar[0, :, :] = np.flipud(list(zip(*grdlat))[::-1])

    grdlonvar = ds.createVariable("merge_grdlon", 'f4', ('time', 'lat', 'lon', ))
    grdlonvar.units = 'unknown'
    grdlonvar[0, :, :] = np.flipud(list(zip(*grdlon))[::-1])

    min = ds.createVariable("merge_mintau", 'f4', ('time', 'lat', 'lon', ))
    min.units = 'unknown'
    min[0, :, :] = np.flipud(list(zip(*mintau))[::-1])

    max = ds.createVariable("merge_maxtau", 'f4', ('time', 'lat', 'lon', ))
    max.units = 'unknown'
    max[0, :, :] = np.flipud(list(zip(*maxtau))[::-1])

    countvar = ds.createVariable("merge_count", 'f4', ('time', 'lat', 'lon', ))
    countvar.units = 'unknown'
    countvar[0, :, :] = np.flipud(list(zip(*count))[::-1])

    sum = ds.createVariable("merge_sumtau", 'f4', ('time', 'lat', 'lon', ))
    sum.units = 'unknown'
    sum[0, :, :] = np.flipud(list(zip(*sumtau))[::-1])
    
    #close file
    ds.close()
    # output = 1 netcdf with six variables (one for each sensor) 
    # merge = run two loops (lat lon) (xdim, ydim from grid) 
    # merge aod equal to average of six sensors making sure no nan
    # identify which sensors are availabe - sensor ID variable 0 or 1 (unavailable vs available)
    # byte type array, keep n-dimensional for now
    # minimum number of things to merge: in the future

    # future capabilities:
    # minimum number, flagging system, use higher quality data (prefiltered)


# filename - 'gridNet.nc' (this is what it saves to)
# sat_files - list of satellite files to pull data from
# saves all satellite sensor gridded data as separate variables
#def grid_nc_single_statistics(limit, gsize, inlat, inlon, geo_list, phy_list, filelist, filename):
def grid_nc_sensor_statistics(limit, gsize, geo_list, phy_list, filelist, filename, phy_nc=None, phy_hdf=None):
    #start_timer = time.time()

    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    lat_dim = int(1+round(((maxlat-minlat)/gsize)))
    lon_dim = int(1+round(((maxlon-minlon)/gsize)))
    
    #set target netcdf file
    ds = netCDF4.Dataset(filename, 'w',format='NETCDF4')
    
    # Create time, lat, lon dimensions
    timed = ds.createDimension('time', None)
    lat = ds.createDimension('lat', lat_dim) #latitude is fatitude (y)
    lon = ds.createDimension('lon', lon_dim) # x
    
    #set time variable to filename parse
    times = ds.createVariable('time', 'f4', ('time',))
    
    #set lats and lons  (or rather y and x)
    lats = ds.createVariable('latitude', 'f4', ('lat',)) #latitude is fatitude
    lons = ds.createVariable('longitude', 'f4', ('lon',))

    #bounds for lats and lons
    lats[:] = np.arange(0,lat_dim)*gsize+minlat
    lons[:] = np.arange(0,lon_dim)*gsize+minlon

    #metadata for lat lon
    
    # unknown number of sensors, save values to list
    values = []
    sensors = {}

    # instantiate empty arrays for summary statistics
    indata_temp = []
    inlat_temp = []
    inlon_temp = []
    
    #run through and assign gridded data to sensor variable
    # filenames are in the form:
    # [[time0: [sat1], [sat2], ... ], ...[timek: [sat1], [sat2], ....]] 
    # we assume we receive a signle dict for a time interval though

    for i, s_name in enumerate(filelist.keys()): 
        start_timer = time.time()

        for s in filelist.get(s_name):
            print("WORK SENSOR: ", s_name)

            L2FID,GeoID,PhyID = sat_data_input.open_file(str(s))

            # check for hdf files 
            if (str(s).split(".")[-1] == "hdf"):
                geo_list = ['Latitude', 'Longitude']
            else:
                geo_list = ['latitude', 'longitude']

            lat,lon,phy_vars, meta = filter_data.filter_data(GeoID, PhyID, geo_list, phy_list, phy_nc, phy_hdf)
            
            #close/end the file
            try:#netcdf
                L2FID.close()
            except:#hdf
                L2FID.end()

            #print("latitude:" , lat)
            #print("Physical variables: ", phy_vars)

            for count, p_var in enumerate(phy_list):
                # statistics appending for this sensor
                if len(indata_temp) < count+1:
                    indata_temp.append([])
                    inlat_temp.append([])
                    inlon_temp.append([])
                    #print("sdf",indata_temp)
                
                indata_temp[count] = ma.concatenate([indata_temp[count], phy_vars[count]])
                inlat_temp[count] = ma.concatenate([inlat_temp[count], lat[count]])
                inlon_temp[count] = ma.concatenate([inlon_temp[count], lon[count]])

        #grid parameters
        for j, p_vars in enumerate(phy_list):

            #print("Indata: ", indata_temp[j])
            avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,float(gsize),indata_temp[j],inlat_temp[j],inlon_temp[j])

            rotated = list(zip(*avgtau))[::-1]

            index = i * len(phy_list) + j
            name = str(s_name+"_"+p_vars)
            values.append(ds.createVariable(name, 'f4', ('time', 'lat', 'lon', ), fill_value=meta[j]["_FillValue"]))
            
            # metadata
            values[index].units = meta[j]["units"]
            values[index].valid_range = meta[j]["valid_range"]
            #values[index]._FillValue = meta[j]["_FillValue"]
            values[index].long_name = meta[j]["long_name"]
            values[index].scale_factor = meta[j]["scale_factor"]
            values[index].add_offset = meta[j]["add_offset"]
            values[index].Parameter_Type = meta[j]["Parameter_Type"]

            #make sure key is in dict
            if not p_vars in sensors:
                sensors[p_vars] = []
            
            values[index][0, :, :] = np.flipud(rotated) # make it 2dimensional
            check = values[index][0, :, :].flatten()
            print("Final: ", check)
            sensors[p_vars].append(np.flipud(rotated))

            end_timer = time.time()
            print("Sensor: ",s_name ," time: ", end_timer - start_timer)
        
        # reset
        indata_temp = []
        inlat_temp = []
        inlon_temp = []

    # put in gridded statistics
    # avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = gridding.grid(limit,float(gsize),indata_stat,inlat_stat,inlon_stat)
    sensor_id = []
    shape = (len(np.arange(0,lat_dim)*gsize+minlat), len(np.arange(0,lon_dim)*gsize+minlon))
    avg = []

    for k, p_vars in enumerate(phy_list):
        # averages
        try:
            avgtau = np.nanmean(np.array(sensors.get(p_vars)), axis=0 )
            # other statistics nan

            # sensor id x
            sensor_id_x = np.zeros(shape)
            for i in range(len(sensors[p_vars])):
                sensor_id_x = sensor_id_x + (np.array(sensors.get(p_vars)[i])*0+1)
                #print(np.nan_to_num(np.array(sensors.get(p_vars)[i])*0+1, -1))

            avg_name = str("merge_avgtau_" + p_vars)
            avg.append(ds.createVariable(avg_name, 'f4', ('time', 'lat', 'lon', )))
            avg[k][0, :, :] = avgtau

            # metadata
            # naming convention, units should be consistent
            # read from configuration file
            # follow netcdf conventions

            # avg[k].units = meta["units"]
            # avg[k].valid_range = meta["valid_range"]
            # avg[k]._FillValue = meta["_FillValue"]
            # avg[k].long_name = meta["long_name"]
            # avg[k].scale_factor = meta["scale_factor"]
            # avg[k].add_offset = meta["add_offset"]
            # avg[k].Parameter_Type = meta["Parameter_Type"]

            sensor_id_name = str("sensor_id_x_" + p_vars)
            sensor_id.append(ds.createVariable(sensor_id_name, 'f4', ('time', 'lat', 'lon', )))
            sensor_id[k].units = 'unknown'
            sensor_id[k][0, :, :] = sensor_id_x
        except:
            #print("There are: ", len(sensors.get(p_vars)), " sensors for the var")
            print(p_vars)
            print(sensors.get(p_vars))
        

    #close file
    ds.close()
    # output = 1 netcdf with six variables (one for each sensor) 
    # merge = run two loops (lat lon) (xdim, ydim from grid) 
    # merge aod equal to average of six sensors making sure no nan
    # identify which sensors are availabe - sensor ID variable 0 or 1 (unavailable vs available)
    # byte type array, keep n-dimensional for now
    # minimum number of things to merge: in the future

    # future capabilities:
    # minimum number, flagging system, use higher quality data (prefiltered)
    