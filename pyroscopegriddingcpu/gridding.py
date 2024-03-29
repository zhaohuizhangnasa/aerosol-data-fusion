# Gridding

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

# L2 data documentation: 
# - https://oceancolor.gsfc.nasa.gov/docs/format/l2nc/
# - https://oceancolor.gsfc.nasa.gov/docs/format/l2nc_viirs/

import numpy as np
import numpy.ma as ma
from pyroscopegriddingcpu.sat_data_input import *
from pyroscopegriddingcpu.filter_data import *


# limit - [lat1,lat2,lon1,lon2]
# gsize - pixel size
# indata - inputs
# inlat - list of latitudes we map for
# inlon - list of longitudes we map for
def grid(limit,gsize,indata,inlat,inlon): #valid_range
    dx=gsize
    dy=gsize
    
    # define boundary coordinates
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])
    
    dxx = gsize/2.0
        
    # pixel dimensions
    xdim=int(1+round((abs(maxlon-minlon)/dx)))
    ydim=int(1+round((abs(maxlat-minlat)/dy)))
    
    sumtau=np.zeros((xdim,ydim)) # init 2d array map with zeros
    sqrtau=np.zeros((xdim,ydim))
    count=np.zeros((xdim,ydim))
    
    mintau=np.full([xdim,ydim],np.inf) # init 2d array map with defaults
    maxtau=np.full([xdim,ydim],-np.inf) # fill value, do not change min max
    
    avgtau=np.full([xdim,ydim],np.nan) # write nans as fill values rather than nans
    stdtau=np.full([xdim,ydim],np.nan)
    grdlat=np.full([xdim,ydim],np.nan)
    grdlon=np.full([xdim,ydim],np.nan)

    for ii in range(len(indata)):
        #check within bounds
        # indata should be filtered based on range (not 0, 5 but rather valid_range)
        if (inlat[ii]>=minlat-dxx and inlat[ii] <= maxlat+dxx and \
            inlon[ii]>= minlon-dxx and inlon[ii]<= maxlon+dxx): # and indata[ii] >0.0 and indata[ii]<=5.0):

            # pixel coordinates for bounds
            i=int(round((inlon[ii]-minlon)/dx))
            j=int(round((inlat[ii]-minlat)/dy))

            if i >= xdim : i = xdim - 1 
            if j >= ydim : j = ydim - 1
            if i<0: i = 0 
            if j< 0: j = 0

            # do not take nan values into account
            # nans should not be counted
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
            grdlat[i,j]=dy*j+minlat
            if count[i,j] > 0:
                avgtau[i,j]=sumtau[i,j]/count[i,j]
                if count[i,j] == 1:
                    stdtau[i,j]= 0
                if count[i,j] > 1:
                    #para1 = (sqrtau[i,j] - (sumtau[i,j]**2 / count[i,j]))
                    para1 = (sqrtau[i,j] - (avgtau[i,j]**2 * count[i,j]))
                    para1 = para1 / (count[i,j] - 1)
                    if para1 >= 0:
                        stdtau[i,j]=np.sqrt(para1)
                        
                else:
                        sumtau[i,j]=np.nan

    # change none to fill values
    mintau[np.isinf(mintau)] = np.nan
    maxtau[np.isinf(maxtau)] = np.nan

    return avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau

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

        # check for hdf files 
        if (filename.split(".")[-1] == "hdf"):
            geo_list = ['Latitude', 'Longitude']
        else:
            geo_list = ['latitude', 'longitude']

        # valid data points
        lat,lon,phy_vars, metadata = filter_data_nc(GeoID, PhyID, geo_list, phy_list, phy_nc=phy_list, phy_hdf=phy_list)
        
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


# filelist - list of sensor names we want to grid here
# returns gridded statistics
def multi_sensor_grid(filelist, gsize, limit, phy_list, geo_list=['latitude', 'longitude']):
    # get gridding parameters
    inlat,inlon,indata = multi_sensor_grid_data(filelist, phy_list, geo_list)
    
    #return gridded statistics
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,gsize,indata,inlat,inlon)
    
    return avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau

#grid(limit,gsize,indata,inlat,inlon)
if __name__ == '__main__':
    limit = [-90., 90., -180., 180.]
    limit= [1, 10, 1, 10]
    limit = [1, 2, 1, 3]
    gsize = 1# 0.25
    indata = [0] * 10
    inlat = [1,2,3,4,5,6,7,8,9,90]
    inlon = [1,2,3,4,5,6,7,8,9,180]
    
    avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau = grid(limit,gsize,indata,inlat,inlon)
    #print(avgtau)
    
    print(avgtau)
    
