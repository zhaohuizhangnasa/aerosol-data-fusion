
from pysolar.solar import get_altitude
import time
import math

import datetime
import numpy as np
import pandas as pd

"""
def get_SZA(row):
    time_obj = pd.to_datetime(row['Date_UTC']) + pd.to_timedelta(row['Time_UTC'], unit='h')
    time_obj = time_obj.replace(tzinfo=dt.timezone.utc)
    lat, lon = row['Lat'], row['Lon']
    
    return get_altitude(lat, lon, time_obj)
"""



# using pysolar library

# calculates solar zenith for singular point
# uses pysolar library
def get_SZA(lat, lon, time_start, time_diff):
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    
    return get_altitude(lat, lon, time_obj.to_pydatetime())

def get_SZA_array_pysolar(limit, gsize, time_start, time_diff):
    # limit = [min lat, max lat, min lon, max lon]
    max_lat = int(1+round(((limit[1]-limit[0])/gsize)))
    max_lon = int(1+round(((limit[3]-limit[2])/gsize)))
    solar_zenith = np.zeros((max_lat, max_lon))
    
    #loop through array
    for i in range(0, max_lat):
        for j in range(0, max_lon):
            lat = limit[0] +  i * gsize
            lon = limit[2] + j * gsize
            
            sza = get_SZA(lat, lon, time_start, time_diff)
            solar_zenith[i][j] = sza

    return solar_zenith


def get_SZA_array(limit, gsize, time_start, time_diff):
    # limit = [min lat, max lat, min lon, max lon]
    max_lat = int(1+round(((limit[1]-limit[0])/gsize)))
    max_lon = int(1+round(((limit[3]-limit[2])/gsize)))
            
    solar_zenith = [[get_SZA(limit[0] +  i * gsize, limit[2] + j * gsize, time_start, time_diff) for j in range(0, max_lon)] for i in range(0, max_lat)]

    return solar_zenith

def lat_lon_array(limit, gsize):
    max_lat = int(1+round(((limit[1]-limit[0])/gsize)))
    max_lon = int(1+round(((limit[3]-limit[2])/gsize)))
    arr = [[[limit[0] +  i * gsize, limit[2] + j * gsize] for j in range(0, max_lon)] for i in range(0, max_lat)]
    
    return arr

def get_altitude_tup(lat,lon, when):
    return get_altitude(lat,lon, when)

def get_SZA_vec(limit, gsize, time_start, time_diff):
    #3d arr: 2d arr of [lat, lon]
    max_lat = int(1+round(((limit[1]-limit[0])/gsize)))#int(1+round(((limit[1]-limit[0])/gsize)))
    max_lon = int(1+round(((limit[3]-limit[2])/gsize)))#int(1+round(((limit[3]-limit[2])/gsize)))
    lats = [[limit[0] +  i * gsize for j in range(0, max_lon)] for i in range(0, max_lat)]
    lons = [[limit[2] +  j * gsize for j in range(0, max_lon)] for i in range(0, max_lat)]
    
    # time object: when
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    time_obj = time_obj.to_pydatetime()
    
    get_alt_vec = np.vectorize(get_altitude_tup, excluded=['when'])
    
    return get_alt_vec(lat=lats, lon=lons, when=time_obj)
    
    

if __name__ == '__main__':
    limit = [-90, 90, -180, 180]
    #limit = [-10, 10, -20, 20]
    gsize = 0.25
    time_start = "2020-01-01 00:00:00"
    time_diff = 30
    
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    print(time_obj)
    
    # pysolar library
    print("\nPYSOLAR list comprehension:\n")
    start_time = time.time()
    sza = get_SZA_array(limit, gsize, time_start, time_diff)
    #print(sza)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    
    print("\nSZAvectorized:\n")
    start_time = time.time()
    a = get_SZA_vec(limit, gsize, time_start, time_diff)
    #print(a)
    print("--- %s seconds ---" % (time.time() - start_time))
