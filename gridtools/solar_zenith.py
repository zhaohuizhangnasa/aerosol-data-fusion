
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

#global variables 
rad = 180/math.pi

### solar declination angle
"""

theta_0 = double(360.0 * double(D - 1) / 365.0) / rad

declin  = 0.396372 - 22.91327*COS(theta_0) 			$
	  + 4.02543*SIN(theta_0) 				$
	  - 0.387205*COS(2.0*theta_0) 				$
	  + 0.051967*SIN(2.0*theta_0) 				$
	  - 0.154527*COS(3.0*theta_0) 				$
	  + 0.084798*SIN(3.0*theta_0)

declin = declin / rad	;convert to radians
"""
def solar_declination_angle(D):
    theta_0 = (360.0 * (D - 1) / 365.0) / rad #converts to radians

    # theta_0 is currently in radians
    declin  = (0.396372 - 22.91327*math.degrees(math.cos(theta_0))
            + 4.02543*math.degrees(math.sin(theta_0) )
            - 0.387205*math.degrees(math.cos(2.0*theta_0))
            + 0.051967*math.degrees(math.sin(2.0*theta_0) )
            - 0.154527*math.degrees(math.cos(3.0*theta_0) )
            + 0.084798*math.degrees(math.sin(3.0*theta_0)))
    
    declin  = (0.396372 - 22.91327*math.cos(theta_0)
            + 4.02543*math.sin(theta_0) 
            - 0.387205*math.cos(2.0*theta_0)
            + 0.051967*math.sin(2.0*theta_0) 
            - 0.154527*math.cos(3.0*theta_0) 
            + 0.084798*math.sin(3.0*theta_0))
    declin = declin / rad  #converts to radians
    return theta_0, declin
    
### Time correction for solar angle
"""    
correct  = 0.004297 + 0.107029*COS(theta_0) 			$
    - 1.837877*SIN(theta_0) - 0.837378*COS(2.0*theta_0) 	$
    - 2.342824*SIN(2.0*theta_0)

angle = (time - 12.0) * 15.0 + lon + correct
IF angle GT  180.0 THEN angle = angle - 360.0
IF angle LT -180.0 THEN angle = angle + 360.0

lat_rad = lat / rad
lon_rad = lon / rad

angle_rad  = angle / rad
"""
def time_correction_solar_angle(theta_0, time, lat, lon): #assume theta_0 in radians
    
    correct  = (0.004297 + 0.107029*math.degrees(math.cos(theta_0)) 
                - 1.837877*math.degrees(math.sin(theta_0))
                - 0.837378*math.degrees(math.cos(2.0*theta_0))
                - 2.342824*math.degrees(math.sin(2.0*theta_0)))
    
    correct = (0.004297 + 0.107029*math.cos(theta_0)
                - 1.837877*math.sin(theta_0)
                - 0.837378*math.cos(2.0*theta_0)
                - 2.342824*math.sin(2.0*theta_0))
                    
                
    angle = (time - 12.0) * 15.0 + lon + correct
    
    if angle > 180:
        angle = angle - 360
    if angle < -180:
        angle = angle + 360
        
    lat_rad = lat / rad
    lon_rad = lon / rad

    angle_rad  = angle / rad
    
    return angle, angle_rad, lat_rad, lon_rad

### Sun zenith
"""
tmp1 = SIN(lat_rad)*SIN(declin) + COS(lat_rad)*COS(declin)*COS(angle_rad)
tmp2   = ABS(tmp1)

IF tmp2 GT 1.1 THEN BEGIN
    PRINT, 'Error in acos argument in sun zenith'
    PRINT, tmp1
ENDIF ELSE BEGIN							$
    IF tmp2 GT 1.0 THEN BEGIN
	IF tmp1 GT 0.0 THEN tmp1=1.0
	IF tmp1 LT 0.0 THEN tmp1=-1.0
    ENDIF
ENDELSE

theta_rad = ACOS(tmp1)
"""
def sun_zenith(lat_rad, lon_rad, angle_rad, declin):
    
    #tmps are in radians
    tmp1 = math.sin(lat_rad)*math.sin(declin) + math.cos(lat_rad)*math.cos(declin)*math.cos(angle_rad)
    #tmp1 = (math.degrees(math.sin(lat_rad))* math.degrees(math.sin(declin)) 
    #        + math.degrees(math.cos(lat_rad)) * math.degrees(math.cos(declin)) * math.degrees(math.cos(angle_rad)))
    tmp2   = abs(tmp1)
    
    if tmp2 > 1.1:
        print("Error in acos argument in sun zenith: ", tmp1)
    else:
        if tmp2 > 1.0:
            if tmp1 > 0:
                tmp1 = 1.0
            if tmp1 < 0:
                tmp1 = -1.0
                
    theta_rad = math.acos(tmp1)
    
    return theta_rad

### Sun Azumith
"""
tmp1 = SIN(ABS(angle_rad)) *COS(declin) / SIN(theta_rad)
azimuth_rad  = ASIN(tmp1)

IF lat GT declin THEN azimuth_rad = 180.0 / rad - azimuth_rad
IF angle  GT 0.0  THEN azimuth_rad = 360.0 / rad - azimuth_rad
"""
def sun_azimuth(lat, lon, angle, angle_rad, theta_rad, declin):
    tmp1 = math.sin(abs(angle_rad)) * math.cos(declin) / math.sin(theta_rad)
    azimuth_rad  = math.asin(tmp1)
    
    if lat > declin:
        azimuth_rad = 180.0 / rad - azimuth_rad
    if angle > 0.0:
        azimuth_rad = 360.0 / rad - azimuth_rad
        
    return azimuth_rad

#converts datetime object to numbers
def date_to_num(datetime_obj):
    D = datetime_obj.timetuple().tm_yday
    time = datetime_obj.hour + (datetime_obj.minute / 60.0) + (datetime_obj.second / 3600)

    return D, time

# calculates solar zenith for singular point
def get_SZA_new(lat, lon, time_start, time_diff=30):
    datetime_obj = pd.to_datetime(time_start)
    datetime_obj = datetime_obj.replace(tzinfo=datetime.timezone.utc)
    
    D, time = date_to_num(datetime_obj)
    
    theta_0, declin = solar_declination_angle(D)
    angle, angle_rad, lat_rad, lon_rad = time_correction_solar_angle(theta_0, time, lat, lon)
    theta_rad = sun_zenith(lat_rad, lon_rad, angle_rad, declin)
    #azimuth_rad = sun_azimuth(lat, lon, angle, angle_rad, theta_rad, declin)
    
    return math.degrees(theta_rad)

# uses pysolar library
def get_SZA(lat, lon, time_start, time_diff):
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    
    #return get_altitude(lat, lon, time_obj.to_pydatetime())
    return (90-get_altitude(lat, lon, time_obj.to_pydatetime()))

def get_SZA_old_manual(latitude, longitude, time_start, time_diff=30):
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
    
    #testing manual pysolar calculations
    latitude = 39.01#37.7749 # San Francisco latitude
    longitude = -77.01#-122.4194 # San Francisco longitude
    datetime_obj = datetime.datetime(2023, 2, 23, 18, 14, 0) # February 17, 2023 at 12:00 PM

    theta = get_SZA(latitude, longitude, datetime_obj)
    print("Solar zenith angle (manual):", theta, "degrees")
    datetime_obj = datetime_obj.replace(tzinfo=datetime.timezone.utc)
    print("Solar zenith angle (pysolar): ", 90 - get_altitude(latitude, longitude, datetime_obj))
    
    