__author__ = "Sally Zhao"
__copyright__ = "Copyright 2023, Pyroscope"
__credits__ = ["Neil Gutkin", "Jennifer Wei", "Pawan Gupta", "Robert Levy", "Xiaohua Pan", "Zhaohui Zhang"]
__version__ = "1.0.0"
__maintainer__ = "Sally Zhao"
__email__ = "zhaosally0@gmail.com"
__status__ = "Production"

# Satellite Data Import
# 
# The satellite data is in the form of netCDF data and will be assumed to be as such going forward.
# 
# Here are the functions:
# 
# 1. open_file : opens the netCDF4 file, expecting a geolocation and geophysical group
# 2. read_geo_data : reads the open_file's return geo 
# 3. read_phy_data : reads the open_file's return phy 

# imports
import netCDF4
import os
import sys
import numpy as np

import pyhdf
from pyhdf import SD

# Opens file through file path
# Reads in two separate groups, geolocation and geophysical
# Opens file through file path
# Reads in two separate groups, geolocation and geophysical
def open_file(L2FileName):
   if not os.path.isfile(L2FileName):
      err = "Error: file not exis {} ... ".format(L2FileName)
      print(err)

      return(False)
   try: # netcdf
      #print("OPENED AS NETCDF")
      L2FID = netCDF4.Dataset(L2FileName,'r',format='NETCDF4')
   except:
      """
      try: # hdf
         print("OPENED AS HDF")
         #L2FID = SD.SD(L2FileName)
         L2FID = netCDF4.Dataset(L2FileName,'r',format='HDF') 
      except:
         #print("touching")
         os.touch(L2FileName)
      """
      L2FID = SD.SD(L2FileName)
      print("HDF\n:", L2FID)
      #L2FID = netCDF4.Dataset(L2FileName,'r')#,format='HDF') 
      
   ftype = os.path.splitext(L2FileName)[1]
   GeoGroupID = L2FID
   PhyGroupID = L2FID
    
   if ftype != '.hdf':
      GeoGroupID =  L2FID.groups['geolocation_data']
      PhyGroupID =  L2FID.groups['geophysical_data']

   return(L2FID,GeoGroupID,PhyGroupID)

# nc files
# Parses the geolocation group after file is read
# Variables are latitude and longitude by default only
def read_geo_data_nc(L2GeoGroup, var_list=None):
   if var_list is None:
      var_list = ['latitude', 'longitude']

   GeoGroup = L2GeoGroup.variables      
   geo = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[]}   
   
   for var in var_list:
      geo['data'].append(GeoGroup[var])
      geo['scale_factor'].append(GeoGroup[var].scale_factor)
      geo['add_offset'].append(GeoGroup[var].add_offset)
      geo['_FillValue'].append(GeoGroup[var]._FillValue)
         
   return geo

# hdf files
# Parses the geolocation group after file is read
# Variables are latitude and longitude by default only
def read_geo_data_hdf(L2GeoGroup, var_list=None):
   
   if var_list is None:
      var_list = ['latitude', 'longitude']

   geo = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[]}
   var_list=[]
   var_list.extend(['Latitude', 'Longitude'])

   for var in var_list:
      sds = L2GeoGroup.select(var)
      data = sds.get()
      metadata = sds.attributes()
      
      #make same as convention
      geo['data'].append(data)
      geo['scale_factor'].append(metadata["scale_factor"])
      geo['add_offset'].append(metadata["add_offset"])
      geo['_FillValue'].append(metadata["_FillValue"])
   
   return geo
   """
   if var_list is None:
      var_list = ['Latitude', 'Longitude']

   GeoGroup = L2GeoGroup.variables      
   geo = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[]}   
   
   for var in var_list:
      geo['data'].append(GeoGroup[var])
      geo['scale_factor'].append(GeoGroup[var].scale_factor)
      geo['add_offset'].append(GeoGroup[var].add_offset)
      geo['_FillValue'].append(GeoGroup[var]._FillValue)
         
   return geo
   """


"""
# Parses the geolocation group after file is read
# Variables are latitude and longitude by default only
def read_geo_data(L2GeoGroup, var_list=None):
   if var_list is None:
      var_list = ['latitude', 'longitude']

   try: # netcdf file
      GeoGroup = L2GeoGroup.variables      
      geo = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[]}   
      for var in var_list:
         geo['data'].append(GeoGroup[var])
         geo['scale_factor'].append(GeoGroup[var].scale_factor)
         geo['add_offset'].append(GeoGroup[var].add_offset)
         geo['_FillValue'].append(GeoGroup[var]._FillValue)
         
   except: # hdf file
      geo = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[]}
      var_list=[]
      var_list.extend(['Latitude', 'Longitude'])

      for var in var_list:
         sds = L2GeoGroup.select(var)
         data = sds.get()
         metadata = sds.attributes()
         
         #make same as convention
         geo['data'].append(data)
         geo['scale_factor'].append(metadata["scale_factor"])
         geo['add_offset'].append(metadata["add_offset"])
         geo['_FillValue'].append(metadata["_FillValue"])
      
      
   return(geo)
"""
"""
# Parses the geophysical group after file is read
# Variables are aod by default only
def read_phy_data(L2PhyGroup, L2GeoGroup, var_list=None, var_nc=None, var_hdf=None):
   if var_list is None:
      var_list = ['Optical_Depth_Land_And_Ocean']

   #try: # netcdf
   
   phy = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[], "valid_range":[], "long_name":[], "units":[], "Parameter_Type":[]}   
   for i in range(len(var_nc)):
      try:
            PhyGroup = L2PhyGroup.variables
            var = var_nc[i]

            #print("NETCDF VARIABLE: ", var)
            phy['data'].append(PhyGroup[var] )
            #print("Read physical in: ", phy['data'])

            phy['valid_range'].append(PhyGroup[var].valid_range)
            phy['_FillValue'].append(PhyGroup[var]._FillValue)
            phy['long_name'].append(PhyGroup[var].long_name)
            phy['units'].append(PhyGroup[var].units)
            phy['scale_factor'].append(PhyGroup[var].scale_factor)
            phy['add_offset'].append(PhyGroup[var].add_offset)
            phy['Parameter_Type'].append(PhyGroup[var].Parameter_Type)
      except:
            # variable not found in PhyGroup
            #print("Physical Var Error: ", AttributeError)

            try: # look now in geolocation variables
               # mainly for variables such as solar_zenith or scatter_angle
               #print("INDEX?: ", i)
               var = var_nc[i]

               GeoGroup = L2GeoGroup.variables 
               #print("GEO VARIABLE: ", var)

               phy['data'].append(GeoGroup[var] )
               #print("Read in: ", phy['data'])

               phy['valid_range'].append(GeoGroup[var].valid_range)
               phy['_FillValue'].append(GeoGroup[var]._FillValue)
               phy['long_name'].append(GeoGroup[var].long_name)
               phy['units'].append(GeoGroup[var].units)
               phy['scale_factor'].append(GeoGroup[var].scale_factor)
               phy['add_offset'].append(GeoGroup[var].add_offset)
               phy['Parameter_Type'].append(GeoGroup[var].Parameter_Type)

            except: # not here, treat as HDF file
               try: 
                  var = var_hdf[i]
                  
                  sds = L2PhyGroup.select(var)
                  data = sds.get()
                  metadata = sds.attributes()
                  
                  #make same as convention
                  phy['data'].append(data)

                  phy['valid_range'].append(metadata["valid_range"])
                  phy['_FillValue'].append(metadata["_FillValue"])
                  phy['long_name'].append(metadata["long_name"])
                  phy['units'].append(metadata["units"])
                  phy['scale_factor'].append(metadata["scale_factor"])
                  phy['add_offset'].append(metadata["add_offset"])
                  phy['Parameter_Type'].append(metadata["Parameter_Type"])
                  
               except: # varaible not found
                  raise AttributeError("No such variable", var_list[i])      
         
   return(phy)
"""
# nc
# Parses the geophysical group after file is read
# Variables are aod by default only
def read_phy_data_nc(L2PhyGroup, L2GeoGroup, var_list=None, var_nc=None):
   if var_list is None:
      var_list = ['Optical_Depth_Land_And_Ocean']

   phy = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[], "valid_range":[], "long_name":[], "units":[], "Parameter_Type":[]}   
   for i in range(len(var_nc)):
      try:
            PhyGroup = L2PhyGroup.variables
            var = var_nc[i]

            #print("NETCDF VARIABLE: ", var)
            phy['data'].append(PhyGroup[var] )
            #print("Read physical in: ", phy['data'])

            phy['valid_range'].append(PhyGroup[var].valid_range)
            phy['_FillValue'].append(PhyGroup[var]._FillValue)
            phy['long_name'].append(PhyGroup[var].long_name)
            phy['units'].append(PhyGroup[var].units)
            phy['scale_factor'].append(PhyGroup[var].scale_factor)
            phy['add_offset'].append(PhyGroup[var].add_offset)
            phy['Parameter_Type'].append(PhyGroup[var].Parameter_Type)
      except:
            # variable not found in PhyGroup
            #print("Physical Var Error: ", AttributeError)

            try: # look now in geolocation variables
               # mainly for variables such as solar_zenith or scatter_angle
               #print("INDEX?: ", i)
               var = var_nc[i]

               GeoGroup = L2GeoGroup.variables 
               #print("GEO VARIABLE: ", var)

               phy['data'].append(GeoGroup[var] )
               #print("Read in: ", phy['data'])

               phy['valid_range'].append(GeoGroup[var].valid_range)
               phy['_FillValue'].append(GeoGroup[var]._FillValue)
               phy['long_name'].append(GeoGroup[var].long_name)
               phy['units'].append(GeoGroup[var].units)
               phy['scale_factor'].append(GeoGroup[var].scale_factor)
               phy['add_offset'].append(GeoGroup[var].add_offset)
               phy['Parameter_Type'].append(GeoGroup[var].Parameter_Type)
            except:
               print("PHY DATA NOT FOUND: NC")

   return(phy)

# HDF
# Parses the geophysical group after file is read
# Variables are aod by default only
def read_phy_data_hdf(L2PhyGroup, L2GeoGroup, var_list=None, var_hdf=None):
   if var_list is None:
      var_list = ['Optical_Depth_Land_And_Ocean']

   #try: # netcdf
   
   phy = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[], "valid_range":[], "long_name":[], "units":[], "Parameter_Type":[]}   
   
   for i in range(len(var_hdf)):
      try: 
         var = var_hdf[i]
         
         sds = L2PhyGroup.select(var)
         data = sds.get()
         #print("HDF DATA: ", data)
         metadata = sds.attributes()
         
         #make same as convention
         phy['data'].append(data)

         phy['valid_range'].append(metadata["valid_range"])
         phy['_FillValue'].append(metadata["_FillValue"])
         phy['long_name'].append(metadata["long_name"])
         phy['units'].append(metadata["units"])
         phy['scale_factor'].append(metadata["scale_factor"])
         phy['add_offset'].append(metadata["add_offset"])
         phy['Parameter_Type'].append(metadata["Parameter_Type"])
         
         #display data
         print("DATA")
         data2 = np.array(data).flatten()
         print(data2[data2!=-9999])
         
         
         
      except: # varaible not found
         raise AttributeError("No such variable", var_list[i])      
     
   print("HDF:", phy)    
   return(phy)

# nc
# Parses the geophysical group after file is read
# Variables are aod by default only
def read_phy_data_hdf_vnew(L2PhyGroup, L2GeoGroup, var_list=None, var_hdf=None):
   if var_list is None:
      var_list = ['Optical_Depth_Land_And_Ocean']

   phy = {"data":[], "scale_factor":[], "add_offset":[], "_FillValue":[], "valid_range":[], "long_name":[], "units":[], "Parameter_Type":[]}   
   for i in range(len(var_hdf)):
      try:
            PhyGroup = L2PhyGroup.variables
            var = var_hdf[i]

            #print("NETCDF VARIABLE: ", var)
            phy['data'].append(PhyGroup[var] )
            #print("Read physical in: ", phy['data'])

            phy['valid_range'].append(PhyGroup[var].valid_range)
            phy['_FillValue'].append(PhyGroup[var]._FillValue)
            phy['long_name'].append(PhyGroup[var].long_name)
            phy['units'].append(PhyGroup[var].units)
            phy['scale_factor'].append(PhyGroup[var].scale_factor)
            phy['add_offset'].append(PhyGroup[var].add_offset)
            phy['Parameter_Type'].append(PhyGroup[var].Parameter_Type)
      except:
            # variable not found in PhyGroup
            #print("Physical Var Error: ", AttributeError)

            try: # look now in geolocation variables
               # mainly for variables such as solar_zenith or scatter_angle
               #print("INDEX?: ", i)
               var = var_hdf[i]

               GeoGroup = L2GeoGroup.variables 
               #print("GEO VARIABLE: ", var)

               phy['data'].append(GeoGroup[var] )
               #print("Read in: ", phy['data'])

               phy['valid_range'].append(GeoGroup[var].valid_range)
               phy['_FillValue'].append(GeoGroup[var]._FillValue)
               phy['long_name'].append(GeoGroup[var].long_name)
               phy['units'].append(GeoGroup[var].units)
               phy['scale_factor'].append(GeoGroup[var].scale_factor)
               phy['add_offset'].append(GeoGroup[var].add_offset)
               phy['Parameter_Type'].append(GeoGroup[var].Parameter_Type)
            except:
               print("PHY DATA NOT FOUND: NC")

   return(phy)