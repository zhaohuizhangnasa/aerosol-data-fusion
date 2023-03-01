
from pysolar.solar import get_altitude
import time
import math

from sunpy.coordinates import frames
#from sunpy.time import Time
from astropy.coordinates import EarthLocation
from pvlib import solarposition
import pytz


import datetime
import numpy as np
import pandas as pd

import numba as nb
from numba import jit
from numba import cuda

from joblib import Parallel, delayed

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
def get_SZA_old(lat, lon, time_start, time_diff):
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    
    #return get_altitude(lat, lon, time_obj.to_pydatetime())
    return (90-get_altitude(lat, lon, time_obj.to_pydatetime()))

def get_SZA(latitude, longitude, time_start, time_diff=30):
    datetime_obj = pd.to_datetime(time_start)
    #print("TIME:", datetime_obj)
    
    # Calculate the day of the year (n)
    n = datetime_obj.timetuple().tm_yday
    
    #print("days: ", n)

    # Calculate the declination angle of the sun (delta)
    delta = 23.45 * math.degrees(math.sin(math.radians(360 / 365 * (284 + n))))

    # Calculate the difference between the local time and solar time (t)
    t = (4 * (longitude - 15 * datetime_obj.hour) / 60)

    # Calculate the hour angle (H)
    #print("hour", datetime_obj.hour, " minute: ", datetime_obj.minute)
    local_time = (datetime_obj.hour + ((datetime_obj.minute + (datetime_obj.second / 60))/60)) + (longitude/15)
    #print("Local time: ", local_time)
    #print("local time first half: ", (datetime_obj.hour + ((datetime_obj.minute + (datetime_obj.second / 60))/60)))
    B = (n-1) * (360/365)
    E = 229.2*(0.000075 + 0.001868*math.degrees(math.cos(math.radians(B)))- 0.032077*math.degrees(math.sin(math.radians(B))) - 0.014615*math.degrees(math.cos(math.radians(2*B)))- 0.04089 *math.degrees(math.sin(math.radians(2*B))))
    solar_time = (datetime_obj.hour - (longitude / 15.0)) + E
    #print("Solar time", solar_time)
    #H = (15 * (12 - datetime_obj.hour) + t) % 360
    local_solar_time = local_time - solar_time
    H = 15 * (12 -local_solar_time)

    # Calculate the solar zenith angle (theta)
    sin_theta = math.sin(math.radians(latitude)) * math.sin(math.radians(delta)) + math.cos(math.radians(latitude)) * math.cos(math.radians(delta)) * math.cos(math.radians(H))
    theta =90- math.degrees(math.acos(math.radians(sin_theta))) #math.degrees(math.asin(sin_theta))

    return theta/.01


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

#@jit(nopython=True)
def get_SZA_parallelized(limit, gsize, time_start, time_diff,num_cores=1):
    # limit = [min lat, max lat, min lon, max lon]
    max_lat = int(1+round(((limit[1]-limit[0])/gsize)))
    max_lon = int(1+round(((limit[3]-limit[2])/gsize)))
            
    #solar_zenith = Parallel(n_jobs=num_cores)([get_SZA(limit[0] +  i * gsize, limit[2] + j * gsize, time_start, time_diff) for j in range(0, max_lon)] for i in range(0, max_lat))
    #res = Parallel(n_jobs=1)(delayed(math.sqrt)(i**2) for i in range(10))
    res = Parallel(n_jobs = num_cores)( 
            delayed( get_SZA )(limit[0] +  i * gsize, limit[2] + j * gsize, time_start, time_diff) 
                                    for i in range(0, max_lat)
                                    for j in range(0, max_lon)
        )
    
    res = np.reshape(res, (-1, max_lon))
    return res#solar_zenith

if __name__ == '__main__':
    limit = [-90, 90, -180, 180]
    #limit = [-10, 10, -20, 20]
    gsize = 0.25
    time_start = "2020-01-01 00:00:00"
    time_diff = 30
    
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    print(time_obj)
    
    # pysolar library
    """
    print("\nPYSOLAR list comprehension:\n")
    start_time = time.time()
    #sza = get_SZA_array(limit, gsize, time_start, time_diff)
    #print(sza)
    #print("--- %s seconds ---" % (time.time() - start_time))
    
    
    print("\nSZAvectorized:\n")
    start_time = time.time()
    #a = get_SZA_vec(limit, gsize, time_start, time_diff)
    #print(len(a))
    #print(a[1])
    print("--- %s seconds ---" % (time.time() - start_time))
    
    print("\nParallelized")
    start_time = time.time()
    #a = get_SZA_parallelized(limit, gsize, time_start, time_diff, -1)
    print(a)
    print("--- %s seconds ---" % (time.time() - start_time))
    """
    
    #testing manual pysolar calculations
    latitude = 39.01#37.7749 # San Francisco latitude
    longitude = -77.01#-122.4194 # San Francisco longitude
    datetime_obj = datetime.datetime(2023, 2, 23, 18, 14, 0) # February 17, 2023 at 12:00 PM

    theta = get_SZA(latitude, longitude, datetime_obj)
    print("Solar zenith angle (manual):", theta, "degrees")
    datetime_obj = datetime_obj.replace(tzinfo=datetime.timezone.utc)
    print("Solar zenith angle (pysolar): ", get_altitude(latitude, longitude, datetime_obj))
    
    #pvlib
    solpos = solarposition.get_solarposition(datetime_obj, latitude, longitude)
    zenith_angle = 90 - solpos['apparent_zenith']
    print("Solar zenith angle (pvlib):", zenith_angle , "degrees")
    
    """
    #sunpy
    location = EarthLocation(lat=latitude, lon=longitude)

    # Calculate the solar position
    frame = frames.HeliographicStonyhurst(obstime=datetime_obj)
    sun_position = frame.transform_to(frames.Heliocentric(observer=location))

    # Calculate the solar zenith angle
    zenith_angle = 90 - sun_position.transform_to(frames.HeliographicStonyhurst).lon.value
    print("Solar zenith angle (sunpy):", zenith_angle, "degrees")
    """



