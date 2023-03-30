__author__ = "Sally Zhao"
__copyright__ = "Copyright 2023, Pyroscope"
__credits__ = ["Neil Gutkin", "Jennifer Wei", "Pawan Gupta", "Robert Levy", "Xiaohua Pan", "Zhaohui Zhang"]
__version__ = "1.0.0"
__maintainer__ = "Sally Zhao"
__email__ = "zhaosally0@gmail.com"
__status__ = "Production"
# Filter data
# We want to filter out our data to make sure what we're reading is valid and ready to be mapped.
# 
# Luckily, the metadata and attributes of the netcdf4 variables gives us an idea of what we should be looking for through:
# 
# 1. valid_range
# 2. units
# 3. scale_factor

import sat_data_input
import numpy as np

# nc
# given GeoID,PhyID (this is read in with open_file) and their corresponding var_lists (variables we care about)
# output lat, lon, and variables we care about (filtered based on metadata)
# lat lon are n size arrays containing lat, long variables for variable n
# phy_vars is n size array containing netcdf4.variables
def filter_data_nc(GeoID,PhyID,geo_list,phy_list, phy_nc = None):
    if geo_list is None:
       geo_list = ["latitude", "longitude"]
    if phy_list is None:
        phy_list = [] 
    
    #basic read in
    geo = sat_data_input.read_geo_data_nc(GeoID, var_list=geo_list)
    phy = sat_data_input.read_phy_data_nc(PhyID, GeoID, phy_list, phy_nc)
    lat = []
    lon = []
    phy_vars = []
    metadata = []
    

    for p in range(len(phy_list)):
        lat1 = geo['data'][0][:,:]
        lon1 = geo['data'][1][:,:]
        phy_vars1 = phy['data'][p][:,:]
        
        #filter based on valid range
        try: # check in the code
            lat1 = lat1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
            lon1 = lon1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
        except:
            #in this case, the phy_variable in question is not mapped to every lat,lon coordinate
            #print("excepted")
            #we try to map it to where data exists for lat and long
            #lat1 = lat1[np.logical_not(np.isnan(phy_vars1))]
            #lon1 = lon1[np.logical_not(np.isnan(phy_vars1))]
            lat1=[]
            lon1=[]
            
        #print("valid range: ", phy['valid_range'][p][0])
        phy_vars1 = phy_vars1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
        
        lat1 = np.array(lat1) * float(geo['scale_factor'][0]) + float(geo['add_offset'][0])
        lat.append(lat1)
        lon1 = np.array(lon1) * float(geo['scale_factor'][1]) + float(geo['add_offset'][1])
        lon.append(lon1)
        #phy_vars1 = np.array(phy_vars1) * float(phy['scale_factor'][p]) + float(phy['add_offset'][p])
        #phy_vars1[phy_vars1 == 0] = 0.5
        phy_vars.append(phy_vars1)

        #print("phyvars1: ", len(phy_vars1))
        
        #print(phy_vars1)

        # metadata attributes
        m = {}
        for key in phy:
            if key == "scale_factor" or key == "add_offset":
                m[key] = float(phy[key][p])
            else:
                #store
                m[key] = phy[key][p]

        metadata.append(m)

    # scale factors
    phy_vars = np.array(phy_vars) #* float(phy['scale_factor'][0]) + float(phy['add_offset'][0])
    #phy_vars = (np.array(phy_vars) - float(phy['add_offset'][0]) )/float(phy['scale_factor'][0])
    for p in phy_vars:
        pass
        #print("FINAL MAXES NC:", p.max())
    
    #print("latitude:", len(lat))
    #print("Physical:", len(phy_vars))
    print("Metadata: ", metadata)

    # probably should put in a check that lat lon phys has same length
    # sometimes phys var may differ between aerosol optical depth and angle
    # ALL variables must pass the filter test otherwise get rid of them

    return lat,lon,phy_vars, metadata


# hdf
# given GeoID,PhyID (this is read in with open_file) and their corresponding var_lists (variables we care about)
# output lat, lon, and variables we care about (filtered based on metadata)
# lat lon are n size arrays containing lat, long variables for variable n
# phy_vars is n size array containing netcdf4.variables
def filter_data_hdf(GeoID,PhyID,geo_list,phy_list, phy_hdf = None):
    if geo_list is None:
       geo_list = ["latitude", "longitude"]
    if phy_list is None:
        phy_list = [] 
    
    #basic read in
    geo = sat_data_input.read_geo_data_hdf(GeoID, var_list=geo_list)
    phy = sat_data_input.read_phy_data_hdf(PhyID, GeoID, phy_list, phy_hdf)
    lat = []
    lon = []
    phy_vars = []
    metadata = []

    for p in range(len(phy_list)):
        lat1 = geo['data'][0][:,:]
        lon1 = geo['data'][1][:,:]
        phy_vars1 = phy['data'][p][:,:]
        
        #filter based on valid range
        try: # check in the code
            lat1 = lat1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
            lon1 = lon1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
        except:
            #in this case, the phy_variable in question is not mapped to every lat,lon coordinate
            #print("excepted")
            #we try to map it to where data exists for lat and long
            #lat1 = lat1[np.logical_not(np.isnan(phy_vars1))]
            #lon1 = lon1[np.logical_not(np.isnan(phy_vars1))]
            lat1=[]
            lon1=[]
            
        #print("valid range: ", phy['valid_range'][p][0])
        print("VARIABLE:", phy_list[p])
        print("VALID RANGE: ",phy['valid_range'][p][0], " to ", phy['valid_range'][p][1])
        phy_vars1 = phy_vars1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
        
        #print("AFTER FILTER MAX: ", phy_vars1.max())
        #print("AFTER FILTER MIN: ", phy_vars1.min())
        
        lat1 = np.array(lat1) * float(geo['scale_factor'][0]) + float(geo['add_offset'][0])
        lat.append(lat1)
        lon1 = np.array(lon1) * float(geo['scale_factor'][1]) + float(geo['add_offset'][1])
        lon.append(lon1)
        #phy_vars1 = np.array(phy_vars1) * float(phy['scale_factor'][p]) + float(phy['add_offset'][p])
        #phy_vars1[phy_vars1 == 0] = 0.5
        phy_vars.append(phy_vars1)

        #print("phyvars1: ", len(phy_vars1))
        
        #print(phy_vars1)

        # metadata attributes
        m = {}
        for key in phy:
            # conversion to appropriate values
            if key == "scale_factor" or key == "add_offset":
                m[key] = float(phy[key][p])
            else:
                #store
                m[key] = phy[key][p]
            

        metadata.append(m)

    # scale factors
    phy_vars = np.array(phy_vars) * float(phy['scale_factor'][0]) + float(phy['add_offset'][0])
    for p in phy_vars:
        pass
        #print("FINAL MAXES HDF:", p.max())
    #print("latitude:", len(lat))
    #print("Physical:", len(phy_vars))
    print("Metadata: ", metadata)

    # probably should put in a check that lat lon phys has same length
    # sometimes phys var may differ between aerosol optical depth and angle
    # ALL variables must pass the filter test otherwise get rid of them

    return lat,lon,phy_vars, metadata

"""
# given GeoID,PhyID (this is read in with open_file) and their corresponding var_lists (variables we care about)
# output lat, lon, and variables we care about (filtered based on metadata)
# lat lon are n size arrays containing lat, long variables for variable n
# phy_vars is n size array containing netcdf4.variables
def filter_data(GeoID,PhyID,geo_list,phy_list, phy_nc = None, phy_hdf=None):
    if geo_list is None:
       geo_list = ["latitude", "longitude"]
    if phy_list is None:
        phy_list = [] 
    
    #basic read in
    geo = sat_data_input.read_geo_data(GeoID, var_list=geo_list)
    phy = sat_data_input.read_phy_data(PhyID, GeoID, phy_list, phy_nc, phy_hdf)
    lat = []
    lon = []
    phy_vars = []
    metadata = []

    for p in range(len(phy_list)):
        lat1 = geo['data'][0][:,:]
        lon1 = geo['data'][1][:,:]
        phy_vars1 = phy['data'][p][:,:]
        
        #filter based on valid range
        try: # check in the code
            lat1 = lat1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
            lon1 = lon1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
        except:
            #in this case, the phy_variable in question is not mapped to every lat,lon coordinate
            #print("excepted")
            #we try to map it to where data exists for lat and long
            #lat1 = lat1[np.logical_not(np.isnan(phy_vars1))]
            #lon1 = lon1[np.logical_not(np.isnan(phy_vars1))]
            lat1=[]
            lon1=[]
            
        #print("valid range: ", phy['valid_range'][p][0])
        phy_vars1 = phy_vars1[np.logical_and(phy_vars1>=phy['valid_range'][p][0],phy_vars1<=phy['valid_range'][p][1])]
        
        lat1 = np.array(lat1) * float(geo['scale_factor'][0]) + float(geo['add_offset'][0])
        lat.append(lat1)
        lon1 = np.array(lon1) * float(geo['scale_factor'][1]) + float(geo['add_offset'][1])
        lon.append(lon1)
        #phy_vars1 = np.array(phy_vars1) * float(phy['scale_factor'][p]) + float(phy['add_offset'][p])
        #phy_vars1[phy_vars1 == 0] = 0.5
        phy_vars.append(phy_vars1)

        #print("phyvars1: ", len(phy_vars1))
        
        #print(phy_vars1)

        # metadata attributes
        m = {}
        for key in phy:
            m[key] = phy[key][p]

        metadata.append(m)

    # scale factors
    phy_vars = np.array(phy_vars) * float(phy['scale_factor'][0]) + float(phy['add_offset'][0])
    #print("latitude:", len(lat))
    #print("Physical:", len(phy_vars))
    print("Metadata: ", metadata)

    # probably should put in a check that lat lon phys has same length
    # sometimes phys var may differ between aerosol optical depth and angle
    # ALL variables must pass the filter test otherwise get rid of them

    return lat,lon,phy_vars, metadata
"""