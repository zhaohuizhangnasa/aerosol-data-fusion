
__author__ = "Sally Zhao"
__copyright__ = "Copyright 2023, Pyroscope"
__credits__ = ["Neil Gutkin", "Jennifer Wei", "Pawan Gupta", "Robert Levy", "Xiaohua Pan", "Zhaohui Zhang"]
__version__ = "1.0.0"
__maintainer__ = "Sally Zhao"
__email__ = "zhaosally0@gmail.com"
__status__ = "Production"
# Solar Zenith
#
# Calculates solar zenith and azimuth angles. 
# Contains manual code calculations validated against pysolar
#

#from pysolar.solar import get_altitude
import time
import math

#from sunpy.coordinates import frames
#from sunpy.time import Time
#from astropy.coordinates import EarthLocation
#from pvlib import solarposition
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
def get_SZA(lat, lon, time_start, time_diff=30):
    datetime_obj = pd.to_datetime(time_start)
    datetime_obj = datetime_obj.replace(tzinfo=datetime.timezone.utc)
    
    D, time = date_to_num(datetime_obj)
    
    theta_0, declin = solar_declination_angle(D)
    angle, angle_rad, lat_rad, lon_rad = time_correction_solar_angle(theta_0, time, lat, lon)
    theta_rad = sun_zenith(lat_rad, lon_rad, angle_rad, declin)
    #azimuth_rad = sun_azimuth(lat, lon, angle, angle_rad, theta_rad, declin)
    
    return math.degrees(theta_rad)
"""
# uses pysolar library
def get_SZA(lat, lon, time_start, time_diff):
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    
    #return get_altitude(lat, lon, time_obj.to_pydatetime())
    return (90-get_altitude(lat, lon, time_obj.to_pydatetime()))
""" 

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
    return res #solar_zenith

if __name__ == '__main__':
    limit = [-90, 90, -180, 180]
    #limit = [-10, 10, -20, 20]
    gsize = 0.25
    time_start = "2020-01-01 00:00:00"
    time_diff = 30
    
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    print(time_obj)
    
    #testing manual pysolar calculations
    latitude = 40.01#37.7749 # San Francisco latitude
    longitude = -77.01#-122.4194 # San Francisco longitude
    datetime_obj = datetime.datetime(2023, 2, 1, 10, 14, 0) # February 17, 2023 at 12:00 PM

    theta = get_SZA(latitude, longitude, datetime_obj)
    print("Solar zenith angle (manual):", theta, "degrees")
    datetime_obj = datetime_obj.replace(tzinfo=datetime.timezone.utc)
    #print("Solar zenith angle (pysolar): ", 90 - get_altitude(latitude, longitude, datetime_obj))
    
    