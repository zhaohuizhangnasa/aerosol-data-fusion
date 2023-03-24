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
from numba import cuda
import math
import sat_data_input
import filter_data


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
    
    # pixel dimensions
    xdim=int(1+round((abs(maxlon-minlon)/dx)))
    ydim=int(1+round((abs(maxlat-minlat)/dy)))
    
    sumtau=np.zeros((xdim,ydim)) # init 2d array map with zeros
    sqrtau=np.zeros((xdim,ydim))
    count=np.zeros((xdim,ydim))
    mintau=np.full([xdim,ydim],5000.0) # init 2d array map with defaults
    maxtau=np.full([xdim,ydim],-100.0) # fill value, do not change min max
    avgtau=np.full([xdim,ydim],-9999.0) # write nans as fill values rather than nans
    stdtau=np.full([xdim,ydim],-9999.0)
    grdlat=np.full([xdim,ydim],-9999.0)
    grdlon=np.full([xdim,ydim],-9999.0)

    # Run on CPU if no GPU is available
    if ((len(cuda.list_devices()) == 0 or len(indata) == 0)):
        for ii in range(len(indata)):
            #check within bounds
            # indata should be filtered based on range (not 0, 5 but rather valid_range)
            if (inlat[ii]>=minlat and inlat[ii] <= maxlat and inlon[ii]>= minlon and inlon[ii]<= maxlon): # and indata[ii] >0.0 and indata[ii]<=5.0):
                # pixel coordinates for bounds
                i=int(round((inlon[ii]-minlon)/dx))
                j=int(round((inlat[ii]-minlat)/dy))

                # 0 check
                # disregard if 0
                if indata[ii] != 0:
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
                grdlat[i,j]=dx*j+minlat
                if count[i,j] > 0:
                    avgtau[i,j]=sumtau[i,j]/count[i,j]
                    para1=(1/count[i,j])*(sqrtau[i,j])+(count[i,j])*avgtau[i,j]-2*(avgtau[i,j])*(sumtau[i,j])
                    if para1 > 0:
                        stdtau[i,j]=np.sqrt(para1)
                        
        # change none to fill values
        mintau[mintau==5000]=None
        maxtau[maxtau==-1]=None
        avgtau[avgtau==-1]=None

        return avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau

    # Run on GPU if GPU is available
    # Create device arrays
    d_sumtau = cuda.to_device(sumtau.flatten())
    d_sqrtau = cuda.to_device(sqrtau.flatten())
    d_count = cuda.to_device(count.flatten())
    d_mintau = cuda.to_device(mintau.flatten())
    d_maxtau = cuda.to_device(maxtau.flatten())
    d_avgtau = cuda.to_device(avgtau.flatten())
    d_stdtau = cuda.to_device(stdtau.flatten())
    d_grdlat = cuda.to_device(grdlat.flatten())
    d_grdlon = cuda.to_device(grdlon.flatten())
    d_indata = cuda.to_device(indata)
    d_limit = cuda.to_device(limit)
    d_inlat = cuda.to_device(inlat)
    d_inlon = cuda.to_device(inlon)

    # Run first kernel
    tpb = (256)
    bpg = (int(np.ceil(len(indata) / tpb)))
    grid_gpu_1[bpg, tpb](gsize, ydim, len(indata),\
        d_indata, d_limit, d_inlat, d_inlon, d_sumtau, d_sqrtau, d_count, d_mintau, d_maxtau)

    # Run second kernel
    # Threads per block
    tpb = (16, 16)
    # Blocks per grid
    bpg_x = int(np.ceil(ydim / tpb[0]))
    bpg_y = int(np.ceil(xdim / tpb[1]))
    bpg = (bpg_x, bpg_y)
    grid_gpu_2[bpg, tpb](gsize, minlon, minlat, xdim, ydim,\
        d_grdlon, d_grdlat, d_count, d_avgtau, d_sumtau, d_sqrtau, d_stdtau, d_mintau, d_maxtau)

    # Transfer data to host
    avgtau = d_avgtau.copy_to_host().reshape([xdim,ydim])
    stdtau = d_stdtau.copy_to_host().reshape([xdim,ydim])
    grdlat = d_grdlat.copy_to_host().reshape([xdim,ydim])
    grdlon = d_grdlon.copy_to_host().reshape([xdim,ydim])
    mintau = d_mintau.copy_to_host().reshape([xdim,ydim])
    maxtau = d_maxtau.copy_to_host().reshape([xdim,ydim])
    count = d_count.copy_to_host().reshape([xdim,ydim])
    sumtau = d_sumtau.copy_to_host().reshape([xdim,ydim])

    mintau[mintau==5000]=None
    maxtau[maxtau==-1]=None
    avgtau[avgtau==-1]=None

    return avgtau,stdtau,grdlat,grdlon,mintau,maxtau,count,sumtau

@cuda.jit
def grid_gpu_1(gsize, ydim, indata_length,\
    indata, limit, inlat, inlon, sumtau, sqrtau, count, mintau, maxtau):
    # Establish thread and block indices
    t = cuda.threadIdx.x
    #-------------------------
    d = cuda.blockDim.x
    #-------------------------
    b = cuda.blockIdx.x

    # Compute global array index
    idx = t + b*d

    # Unpack lat/lon limits
    minlat=float(limit[0])
    maxlat=float(limit[1])
    minlon=float(limit[2])
    maxlon=float(limit[3])    

    # Grid values
    if idx < indata_length:
        if (inlat[idx]>=minlat and inlat[idx] <= maxlat and inlon[idx]>= minlon and inlon[idx]<= maxlon): # and indata[ii] >0.0 and indata[ii]<=5.0):
            # pixel coordinates for bounds
            i=int(round((inlon[idx]-minlon)/gsize))
            j=int(round((inlat[idx]-minlat)/gsize))
            ij = j + i*ydim
            
            # do not take nan values into account
            # nans should not be counted
            if indata[idx] != 0:
                cuda.atomic.add(sumtau, ij, indata[idx])
                cuda.atomic.add(sqrtau, ij, indata[idx]**2)
                cuda.atomic.add(count, ij, 1)
                cuda.atomic.min(mintau, ij, indata[idx])
                cuda.atomic.max(maxtau, ij, indata[idx])


@cuda.jit
def grid_gpu_2(gsize, minlon, minlat, xdim, ydim,\
    grdlon, grdlat, count, avgtau, sumtau, sqrtau, stdtau, mintau, maxtau):
    # Establish thread and block indices
    tx = cuda.threadIdx.x
    ty = cuda.threadIdx.y
    #-------------------------
    dx = cuda.blockDim.x
    dy = cuda.blockDim.y
    #-------------------------
    bx = cuda.blockIdx.x
    by = cuda.blockIdx.y

    # Compute global array index
    row = ty + by*dy
    col = tx + bx*dx
    idx = col + row*ydim

    # Grid values
    if row < xdim and col < ydim:
        grdlon[idx] = gsize*row + minlon
        grdlat[idx] = gsize*col + minlat
        if count[idx] > 0:
            avgtau[idx] = sumtau[idx] / float(count[idx])
            para1=(float(1)/float(count[idx]))*(sqrtau[idx])+(count[idx])*avgtau[idx]-2*(avgtau[idx])*(sumtau[idx])
            if para1 > 0:
                stdtau[idx] = math.sqrt(para1)

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
        L2FID,GeoID,PhyID = sat_data_input.open_file(filename)

        # check for hdf files 
        if (filename.split(".")[-1] == "hdf"):
            geo_list = ['Latitude', 'Longitude']
        else:
            geo_list = ['latitude', 'longitude']

        # valid data points
        lat,lon,phy_vars, metadata = filter_data.filter_data(GeoID, PhyID, geo_list, phy_list, phy_nc=phy_list, phy_hdf=phy_list)
        
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
    