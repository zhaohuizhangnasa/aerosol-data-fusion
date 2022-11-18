# %% [markdown]
# # Data fusion - Merging Sensor Data on Variables

# %% [markdown]
# ### Imports

# %%
import pandas as pd
import netCDF4
#import nctoolkit as nc
from netCDF4 import Dataset

import math
import numpy as np 
import numpy.ma as ma
import datetime
import os

import matplotlib
import matplotlib.cm as cm
from matplotlib import pyplot as plt

#mapping
import xarray as xr
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#from cartopy.feature.nightshade import Nightshade
#import shapefile as shp
#import geopandas as gpd

#animation
import matplotlib.animation as animation

#griddign geoTIFF
#from osgeo import gdal
#from osgeo import osr

# %% [markdown]
# ## Satellite Data Import
# 
# Here we use Zhaohui's provided code to read in satellite and sensor data.
# The satellite data is in the form of netCDF data and will be assumed to be as such going forward.
# 
# Here are the functions:
# 
# 1. open_file : opens the netCDF4 file, expecting a geolocation and geophysical group
# 2. read_geo_data : reads the open_file's return geo 
# 3. read_phy_data : reads the open_file's return phy 

# %%
# Opens file through file path
# Reads in two separate groups, geolocation and geophysical
def open_file(L2FileName):
    if not os.path.isfile(L2FileName):
       err = "Error: file not exis {} ... ".format(L2FileName)
       print(err)
       return(False)
    try:
       L2FID = netCDF4.Dataset(L2FileName,'r',format='NETCDF4')
    except IOError:
       os.touch(L2FileName)

    ftype = os.path.splitext(L2FileName)[1]
    GeoGroupID = L2FID
    PhyGroupID = L2FID
    if ftype != '.hdf':
        print(L2FileName)
        GeoGroupID =  L2FID.groups['geolocation_data']
        PhyGroupID =  L2FID.groups['geophysical_data']

    return(L2FID,GeoGroupID,PhyGroupID)

# %%
# Parses the geolocation group after file is read
# Variables are latitude and longitude by default only
def read_geo_data(L2GeoGroup, var_list=None):
    if var_list is None:
       var_list = ['latitude', 'longitude']
    GeoGroup = L2GeoGroup.variables      
    geo = {"data":[], "scalefactor":[], "offset":[], "fillvalue":[]}   
    for var in var_list:
       geo['data'].append(GeoGroup[var])
       geo['scalefactor'].append(GeoGroup[var].scale_factor)
       geo['offset'].append(GeoGroup[var].add_offset)
       geo['fillvalue'].append(GeoGroup[var]._FillValue)
       
    return(geo)

# %%
# Parses the geophysical group after file is read
# Variables are aod by default only
def read_phy_data(L2PhyGroup, var_list=None):
    if var_list is None:
       var_list = ['Optical_Depth_Land_And_Ocean']
 
    PhyGroup = L2PhyGroup.variables
    phy = {"data":[], "scalefactor":[], "offset":[], "fillvalue":[]}   
    for var in var_list:
       phy['data'].append(PhyGroup[var] )
       phy['scalefactor'].append(PhyGroup[var].scale_factor)
       phy['offset'].append(PhyGroup[var].add_offset)
       phy['fillvalue'].append(PhyGroup[var]._FillValue)

    return(phy)

# %% [markdown]
# Reading in files now. Code is Zhaohui's.
# 
# If PRINT_DATA == True, then aod at particular latitude and longitude would be printed out
# 
# I am setting script to false for now as we do not need the contents of this particular cell.

# %% [markdown]
# #### File locations

# %%
# sample file names
fabi16 = "AERDT_L2_ABI_G16.A2021008.1400.001.nc"
fabi17 = "AERDT_L2_ABI_G17.A2021008.1800.001.nc"
fahi08 = "AERDT_L2_AHI_H08.A2021008.0400.001.nc"
fviirs = "AERDT_L2_VIIRS_SNPP.A2021008.2018.001.2021009072852.nc"
fmodis = "MOD04_L2.A2021008.2225.061.2021009072555.hdf"

fpath = "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/" 

# %%
#%%script false 

if __name__ == "__main__":
  filename = fpath+fabi17
  MODIS = False
  print('filename:',filename)

  geo_list = ['latitude', 'longitude']
  phy_list = ['Optical_Depth_Land_And_Ocean']
  if MODIS:
    geo_list = ['Latitude', 'Longitude']
  
  L2FID,GeoID,PhyID = open_file(filename)
  geo = read_geo_data(GeoID, var_list=geo_list)
  phy = read_phy_data(PhyID, var_list=phy_list)
  lat = geo['data'][0][:,:]
  lon = geo['data'][1][:,:]
  aod = phy['data'][0][:,:]

  row,col=np.array(aod).shape
  print("dims:", row, col)

  # valid data points
  lat1=lat[aod>0]
  lon1=lon[aod>0]
  aod1=aod[aod>0]

  #PRINTS
  PRINT_DATA = True
  if PRINT_DATA:
    for i in range(10): #(aod1.size):
      print('(lat,lon)=({0:.2f}, {1:.2f}) aod= {2:.3f}'.\
             format(lat1[i], lon1[i], aod1[i]))
  #print(len(lon1))
  L2FID.close()

# %% [markdown]
# Now we try with variables different from aod

# %%
filename = fpath+fahi08
MODIS = False
print('filename:',filename)

geo_list = ['latitude', 'longitude']
phy_list = ['Optical_Depth_Land_And_Ocean', 'Optical_Depth_Large_Average_Ocean']
if MODIS:
    geo_list = ['Latitude', 'Longitude']
  
L2FID,GeoID,PhyID = open_file(filename)
geo = read_geo_data(GeoID, var_list=geo_list)
phy = read_phy_data(PhyID, var_list=phy_list)
lat = geo['data'][0][:,:]
lon = geo['data'][1][:,:]
aod = phy['data'][0][:,:]

row,col=np.array(aod).shape
print("dims:", row, col)

# valid data points
lat1=lat[aod>0]
lon1=lon[aod>0]
aod1=aod[aod>0]

phy['data'][0].valid_range

aod1
L2FID.close()

# %% [markdown]
# We want to filter out our data to make sure what we're reading is valid and ready to be mapped.
# 
# Luckily, the metadata and attributes of the netcdf4 variables gives us an idea of what we should be looking for through:
# 
# 1. valid_range
# 2. units
# 3. scale_factor

# %%
# given GeoID,PhyID (this is read in with open_file) and their corresponding var_lists (variables we care about)
# output lat, lon, and variables we care about (filtered based on metadata)
# lat lon are n size arrays containing lat, long variables for variable n
# phy_vars is n size array containing netcdf4.variables
def filter_data(GeoID,PhyID,geo_list,phy_list):
    if geo_list is None:
       geo_list = ["latitude", "longitude"]
    if phy_list is None:
        phy_list = []
    
    #basic read in
    geo = read_geo_data(GeoID, var_list=geo_list)
    phy = read_phy_data(PhyID, var_list=phy_list)
    lat = []
    lon = []
    phy_vars = []
    for p in range(len(phy_list)):
        #print("iteration: ", p)
        lat1 = geo['data'][0][:,:]
        lon1 = geo['data'][1][:,:]
        phy_vars1 = phy['data'][p][:,:]
        
        #filter based on valid range
        try:
            lat1 = lat1[np.logical_and(phy_vars1>=phy['data'][p].valid_range[0],phy_vars1<=phy['data'][p].valid_range[1])]
            lon1 = lon1[np.logical_and(phy_vars1>=phy['data'][p].valid_range[0],phy_vars1<=phy['data'][p].valid_range[1])]
        except:
            #in this case, the phy_variable in question is not mapped to every lat,lon coordinate
            #print("excepted")
            #we try to map it to where data exists for lat and long
            #lat1 = lat1[np.logical_not(np.isnan(phy_vars1))]
            #lon1 = lon1[np.logical_not(np.isnan(phy_vars1))]
            lat1=[]
            lon1=[]
            
            
        phy_vars1 = phy_vars1[np.logical_and(phy_vars1>=phy['data'][p].valid_range[0],phy_vars1<=phy['data'][p].valid_range[1])]
        
        lat.append(lat1)
        lon.append(lon1)
        phy_vars.append(phy_vars1)
        

    return lat,lon,phy_vars

# %%
filename = fpath+fahi08
MODIS = False
print('filename:',filename)

geo_list = ['latitude', 'longitude']
phy_list = ['Optical_Depth_Land_And_Ocean', 'Optical_Depth_Large_Average_Ocean', 'Optical_Depth_Ratio_Small_Ocean_0p55micron']
if MODIS:
    geo_list = ['Latitude', 'Longitude']
  
L2FID,GeoID,PhyID = open_file(filename)

row,col=np.array(aod).shape
print("dims:", row, col)

# valid data points
lat,lon,phy_vars = filter_data(GeoID, PhyID, geo_list, phy_list)
#print(len(lat[0]))
#print(phy_vars)

L2FID.close()

# %% [markdown]
# ## Gridding

# %% [markdown]
# To transform the level 2 data into level 3 data, the derived geophysical variables from the L2 data are:
# 
# 1. mapped to uniform space-time grid scales
# 2. weighted based on data quality (for control, consistency, and completeness)
# 
# The L3 data is a summary of a 30 minute interval of L2 data.
# 
# Output format: GeoTIFF
# 
# Resources:
# - https://earthdata.nasa.gov/collaborate/open-data-services-and-software/data-information-policy/data-levels
# - https://earthdata.nasa.gov/learn/articles/gedi-level3data-available
# 

# %% [markdown]
# L2 data documentation: 
# - https://oceancolor.gsfc.nasa.gov/docs/format/l2nc/
# - https://oceancolor.gsfc.nasa.gov/docs/format/l2nc_viirs/

# %%
# code from grdthedata.py

# limit - [lat1,lat2,lon1,lon2]
# gsize - pixel size
# indata - inputs
# inlat - list of latitudes we map for
# inlon - list of longitudes we map for
def grid(limit,gsize,indata,inlat,inlon):
    dx=gsize
    dy=gsize
    
    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    xdim=int(1+round(((maxlon-minlon)/dx)))
    ydim=int(1+round(((maxlat-minlat)/dy)))
    
    sumtau=np.zeros((xdim,ydim)) # init 2d array map with zeros
    sqrtau=np.zeros((xdim,ydim))
    count=np.zeros((xdim,ydim))
    mintau=np.full([xdim,ydim],10.0) # init 2d array map with defaults
    maxtau=np.full([xdim,ydim],-1.0)
    avgtau=np.full([xdim,ydim],-1.0)
    stdtau=np.full([xdim,ydim],-1.0)
    grdlat=np.full([xdim,ydim],-1.0)
    grdlon=np.full([xdim,ydim],-1.0)
    n=0
    
    for ii in range(len(indata)):
        #check within bounds
        #check 0<indata<=5
        if (inlat[ii]>=minlat and inlat[ii] <= maxlat and inlon[ii]>= minlon and inlon[ii]<= maxlon and indata[ii] >0.0 and indata[ii]<=5.0):
            # pixel coordinates for bounds
            i=int(round((inlon[ii]-minlon)/dx))
            j=int(round((inlat[ii]-minlat)/dy))
            
            sumtau[i,j]=sumtau[i,j]+indata[ii]
            sqrtau[i,j]=sqrtau[i,j]+(indata[ii])**2
            count[i,j]+=1
            if indata[ii] < mintau[i,j]:
                mintau[i,j]=indata[ii]
            if indata[ii] > maxtau[i,j]:
                maxtau[i,j]=indata[ii]
                
            #print("updated: ", sumtau[i,j], mintau[i,j], maxtau[i,j])
    for i in range(xdim):
        for j in range(ydim):
            grdlon[i,j]=dx*i+minlon
            grdlat[i,j]=dx*j+minlat
            if count[i,j] > 0:
                avgtau[i,j]=sumtau[i,j]/count[i,j]
                para1=(1/count[i,j])*(sqrtau[i,j])+(count[i,j])*avgtau[i,j]-2*(avgtau[i,j])*(sumtau[i,j])
                if para1 > 0:
                    stdtau[i,j]=np.sqrt(para1)
    mintau[mintau==10]=None
    maxtau[maxtau==-1]=None
    avgtau[avgtau==-1]=None
    return avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau

# %%
fpath = "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/" 
filename = fpath+fabi16
MODIS = False
print('filename:',filename)

geo_list = ['latitude', 'longitude']
phy_list = ['Optical_Depth_Land_And_Ocean', 'Optical_Depth_Large_Average_Ocean', 'Optical_Depth_Ratio_Small_Ocean_0p55micron']
if MODIS:
    geo_list = ['Latitude', 'Longitude']
  
L2FID,GeoID,PhyID = open_file(filename)


# valid data points
lat,lon,phy_vars = filter_data(GeoID, PhyID, geo_list, phy_list)
#print("Phy:", phy_vars[0])
print("Lat:", lat[0])
#print("Lon:", lon)


L2FID.close()

# %%
#grid parameters
limit = [min(lat[0]), max(lat[0]), min(lon[0]),  max(lon[0])]
gsize = 0.25
indata = phy_vars[0]
inlat=lat[0]
inlon = lon[0]

# %%
avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,gsize,indata,inlat,inlon)

# %% [markdown]
# ## Grid to GeoTIFF

# %%
#filename - 'myGeoTIFF.tif'
def create_geoTIFF(limit,gsize,indata,inlat,inlon, filename):
    dx=gsize
    dy=gsize
    
    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    xdim=int(1+round(((maxlon-minlon)/dx)))
    ydim=int(1+round(((maxlat-minlat)/dy)))
    
    #set geotransform
    geotransform = (minlat, dx, 0, maxlon, 0, -dy)
    
    #create n-band raster file
    #each band is a phys_var (number determined by in_data which we put in)
    dst_ds = gdal.GetDriverByName('GTiff').Create(filename, ydim, xdim, len(indata), gdal.GDT_Float32)
    
    dst_ds.SetGeoTransform(geotransform)    # specify coords
    srs = osr.SpatialReference()            # establish encoding
    srs.ImportFromEPSG(3857)                # WGS84 lat/long
    dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
    
    for i,d in enumerate(indata):
        dst_ds.GetRasterBand(i+1).WriteArray(d)   # write r-band to the raster
        
    dst_ds.FlushCache()                     # write to disk
    
    dst_ds = None
    return gdal.Open(filename).ReadAsArray()
    

# %%
fpath = "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/" 
filename = fpath+fabi17
MODIS = False
print('filename:',filename)

geo_list = ['latitude', 'longitude']
#optical land and ocean - one parameter at a time
phy_list = ['Optical_Depth_Land_And_Ocean', 'Optical_Depth_Large_Average_Ocean', 'Optical_Depth_Ratio_Small_Ocean_0p55micron']
if MODIS:
    geo_list = ['Latitude', 'Longitude']
  
L2FID,GeoID,PhyID = open_file(filename)


# valid data points
lat,lon,phy_vars = filter_data(GeoID, PhyID, geo_list, phy_list)
#print("Phy:", phy_vars[0])
print("Lat:", lat[0])
#print("Lon:", lon)


L2FID.close()

# %%
#grid parameters
limit = [min(lat[0]), max(lat[0]), min(lon[0]),  max(lon[0])]
gsize = 0.25
indata = phy_vars[0]
inlat=lat[0]
inlon = lon[0]

avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,gsize,indata,inlat,inlon)

#sample
d = create_geoTIFF(limit,gsize,[mintau],inlat,inlon, fpath+"geo_test_aerdt_v13.tif")

# %%


# %% [markdown]
# ## Color Grids

# %%
#filename - 'myGeoTIFF.tif' with color
def create_geoTIFF_color(limit,gsize,indata,inlat,inlon, filename):

    dx=gsize
    dy=gsize
    
    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    # pixel dimensions
    xdim=int(1+round(((maxlon-minlon)/dx)))
    ydim=int(1+round(((maxlat-minlat)/dy)))
    
    #set geotransform
    geotransform = (minlat, dx, 0, maxlon, 0, -dy)
    
    #create n-band raster file
    #each band is a phys_var (number determined by in_data which we put in)
    dst_ds = gdal.GetDriverByName('GTiff').Create(filename, ydim, xdim, len(indata), gdal.GDT_Int32)
    
    dst_ds.SetGeoTransform(geotransform)    # specify coords
    srs = osr.SpatialReference()            # establish encoding
    srs.ImportFromEPSG(3857)                # WGS84 lat/long
    dst_ds.SetProjection(srs.ExportToWkt()) # export coords to file
    
    #color
    colors = gdal.ColorTable()
    d = indata[0]
    min_start = int(min(d[d!=-10000]))
    max_end = int(max(d[d!=-10000]))
    print(min_start)
    print(max_end)
    
    colors.CreateColorRamp(-0, (255, 0, 0), 10, (242, 0, 0))

    #colors.CreateColorRamp(-10000, (255, 255, 255), int(min_start)-1, (255, 255, 255))
    #colors.CreateColorRamp(min_start, (255,0, 0), int(max_end-min_start/4) + min_start, (125, 125, 0))
    #colors.CreateColorRamp(int(max_end-min_start/4) + min_start, (125,125, 0), int(max_end-min_start/2) + min_start, (0, 255, 0))
    #colors.CreateColorRamp(int(max_end-min_start/2) + min_start, (0,255, 0), int(max_end-min_start*3/4) + min_start, (0, 125, 125))
    #colors.CreateColorRamp(int(max_end-min_start*3/4) + min_start, (125,125, 0), max_end, (0, 0, 255))
    
    for i,d in enumerate(indata):
        dst_ds.GetRasterBand(i+1).WriteArray(d)   # write r-band to the raster\
        
        band = dst_ds.GetRasterBand(i+1)
    
        band.SetRasterColorTable(colors)
        band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
        del band
        
    dst_ds.FlushCache()                     # write to disk
    
    dst_ds = None
    #return gdal.Open(filename).ReadAsArray()
    

# %%
d = create_geoTIFF_color(limit,gsize,[mintau],inlat,inlon, fpath+"geo_test_aerdt_v9_color.tif")

# %% [markdown]
# ## Multiple Sensor Grid
# 
# We want to grid multiple sensors at once. 
# 
# Here is the logic:
# 
# Gridding takes in: 
# - limit - [lat1,lat2,lon1,lon2] [lon1, lon2, lat1, lat2]
# - gsize - pixel size
# - indata - inputs
# - inlat - list of latitudes we map for
# - inlon - list of longitudes we map for
# 
# To grid multiple sensors at once, limit and gsize remain the same.
# 
# However, we append every sensors indata, inlat, and inlon together.
# 
# The separate individual indata, inlat, and inlon should be the same size (n). Together, we will have three 1D arrays of 3n to grid.

# %%
# filelist - list of sensor names we want to grid here
# modis - true false variables that tells us if a sensor is Modis or not
# geolist - geolocation variables to look for
# phylist - list of physical variables we filter
# returns merged list of data, lon, lat we want to grid
def multi_sensor_grid_data(filelist, phy_list, geo_list=['latitude', 'longitude']):
    # instantiate empty arrays
    indata = []
    inlat = []
    inlon = []
    
    #open files
    #append data in
    for filename in filelist:
        L2FID,GeoID,PhyID = open_file(filename)

        # valid data points
        lat,lon,phy_vars = filter_data(GeoID, PhyID, geo_list, phy_list)
        
        #append
        #indata = np.concatenate(indata, np.ndarray(phy_vars[0]).flatten()).flatten()
        #inlat = np.concatenate(inlat, np.ndarray(lat[0]).flatten()).flatten()
        #inlon = np.concatenate(inlon, np.ndarray(lon[0]).flatten()).flatten()
        
        if len(indata) == 0:
            indata = phy_vars[0]
            inlat = lat[0]
            inlon = lon[0]
        else:
            inlat = ma.concatenate([inlat, lat[0]])
            indata = ma.concatenate([indata, phy_vars[0]]) 
            inlon = ma.concatenate([inlon, lon[0]])

        L2FID.close()
    
    return inlat, inlon, indata

# %%
# filelist - list of sensor names we want to grid here
# returns gridded statistics
def multi_sensor_grid(filelist, gsize, limit, phy_list, geo_list=['latitude', 'longitude']):
    # get gridding parameters
    inlat,inlon,indata = multi_sensor_grid_data(filelist, phy_list, geo_list)
    
    #return gridded statistics
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,gsize,indata,inlat,inlon)
    
    return avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau
    

# %%
# multi grid example
fpath = "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/"
fabi16 = "AERDT_L2_ABI_G16.A2021008.1400.001.nc"
fabi17 = "AERDT_L2_ABI_G17.A2021008.1800.001.nc"
fahi08 = "AERDT_L2_AHI_H08.A2021008.0400.001.nc"
fviirs = "AERDT_L2_VIIRS_SNPP.A2021008.2018.001.2021009072852.nc"

sensors = [#fpath + fabi16,
           fpath + fabi17,
           fpath + fahi08,
           fpath  + fviirs]

# %%
L2FID,GeoID,PhyID = open_file(fpath + fabi17)
L2FID.close()

# %%
geo_list = ['latitude', 'longitude']
phy_list = ['Optical_Depth_Land_And_Ocean']
lat,lon,phy_vars = multi_sensor_grid_data(sensors, phy_list, geo_list)

# %%
# gridding parameters
limit = [min(lat), max(lat), min(lon),  max(lon)]
gsize = 0.25 
indata = phy_vars
inlat=lat
inlon = lon

# %%
avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau =  multi_sensor_grid(sensors, gsize, limit, phy_list, geo_list)

# %%
# output this multi-sensor grid to a geotiff
#sample
d2 = create_geoTIFF(limit,gsize,[mintau],inlat,inlon, fpath+"odla_aerdt_merged_v2.tif")

# %% [markdown]
# ## Misc Helper Functions
# 
# 1. filename_time(filename) 
# - filename does not include path
# - returns datetime object (ex. datetime.datetime(2021, 1, 8, 8, 50))

# %%
# due to the introduction of the time variable, we need to know the time of the sensor data
# this information is not given in the file itself, but rather the file name
# e.g. "AERDT_L2_ABI_G16.A2021008.0850.001.nc" --> 10/08/2021 8:50 AM

# function that pulls time from filename (str)
# returns datetime object
def filename_time(filename):
    time = filename[-17:-7] #following convention
    date = datetime.datetime.strptime(time.split(".")[0], '%y%j')
    day_time = time.split(".")[1]
    date = date.replace(hour=int(day_time[:2]), minute=int(day_time[2:]))
    
    return date

# %%
#Example return
filename_time("AERDT_L2_ABI_G16.A2021008.0850.001.nc")

# %% [markdown]
# ## Saving as NetCDF
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

# %%
# geoTiff to netcdf
# given function
fpath1 = "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/" 
fpath2 = "C:/Users/swzhao/Desktop/Data exports/"
ds = gdal.Translate(fpath2 + "geotiffv13_translate.nc", fpath1+"geo_test_aerdt_v13.tif", format='NetCDF')

# %%
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
    
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,gsize,indata,inlat,inlon)
    
    rotated = list(zip(*mintau))[::-1]
    value[0, :, :] = (np.flipud(rotated))
    
    #close file
    ds.close()
    

# %%
#grid parameters
limit = [min(lat[0]), max(lat[0]), min(lon[0]),  max(lon[0])]
gsize = 0.1
indata = phy_vars[0]
inlat=lat[0]
inlon = lon[0]

# %%
fpath = "C:/Users/swzhao/Desktop/Data exports/"
gridNC_single(limit, gsize, indata, inlat, inlon, fpath+"test_random_grid_19.nc")

# %%
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
    """
    for i in range(0,len(sat_files),time_interval):
        satellites = sat_files[i: i+time_interval]
        print(satellites)
        avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau =  multi_sensor_grid(satellites, gsize, limit, phy_list, geo_list)
    
        rotated = list(zip(*mintau))[::-1]
        value[i, :, :] = (np.flipud(rotated))
        """
    
    """
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau =  multi_sensor_grid(sat_files, gsize, limit, phy_list, geo_list)
    
    rotated = list(zip(*mintau))[::-1]
    value[0, :, :] = (np.flipud(rotated))
    """
    
    for i in range(0,len(sat_files),time_interval):
        satellites = sat_files[i: i+time_interval]
        print(satellites)
        avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau =  multi_sensor_grid(satellites, gsize, limit, phy_list, geo_list)
    
        rotated = list(zip(*mintau))[::-1]
        value[i, :, :] = (np.flipud(rotated))
        
    
    #close file
    ds.close()
    

# %%
#file paths for all in folder
fpath = "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/ABI_G16/"
sat_files = [fpath+f for f in os.listdir(path=fpath)]
sat_files

# %%
# multi grid example
fpath = "C:/Users/swzhao/Desktop/Old Code and Satellite Inputs/Satellite Sensor Data/"
fabi16 = "AERDT_L2_ABI_G16.A2021008.1400.001.nc"
fabi17 = "AERDT_L2_ABI_G17.A2021008.1800.001.nc"
fahi08 = "AERDT_L2_AHI_H08.A2021008.0400.001.nc"
fviirs = "AERDT_L2_VIIRS_SNPP.A2021008.2018.001.2021009072852.nc"

sat_files = [#fpath + fabi16,
           fpath + fabi17,
           fpath + fahi08,
           fpath  + fviirs]

# %%
geo_list = ['latitude', 'longitude']
phy_list = ['Optical_Depth_Land_And_Ocean']
lat,lon,phy_vars = multi_sensor_grid_data(sat_files, phy_list, geo_list)

# %%
# gridding parameters
limit = [min(lat), max(lat), min(lon),  max(lon)]
gsize = 1.0 
indata = phy_vars
inlat=lat
inlon = lon

# %%
fpath = "C:/Users/swzhao/Desktop/Data exports/"
gridNC_time(limit, gsize, indata, inlat, inlon, 1, phy_list, geo_list, sat_files, fpath + "NCtimetestv7.nc")

# %%



