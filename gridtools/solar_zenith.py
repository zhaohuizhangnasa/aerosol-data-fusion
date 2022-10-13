
from pysolar.solar import get_altitude
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


def into_range(x, range_min, range_max):
    shiftedx = x - range_min
    delta = range_max - range_min
    return (((shiftedx % delta) + delta) % delta) + range_min

# Returns azumith, elevation of sun for singular time instance
# when - date variable (year, month, day, hour, minute, second, timezone)
# location - lat, lon tuple (lat, lon)
# refraction - boolean
def sun_position(when, location, refraction):
    # Extract the passed data
    year = when.year
    month = when.month
    day = when.day
    hour = when.hour
    minute = when.minute
    second = when.second
    timezone = 4 #default utc
    #year, month, day, hour, minute, second, timezone = when
    latitude, longitude = location

    # Math typing shortcuts
    rad, deg = math.radians, math.degrees
    sin, cos, tan = math.sin, math.cos, math.tan
    asin, atan2 = math.asin, math.atan2

    # Convert latitude and longitude to radians
    rlat = rad(latitude)
    rlon = rad(longitude)

    # Decimal hour of the day at Greenwich
    greenwichtime = hour - timezone + minute / 60 + second / 3600

    # Days from J2000, accurate from 1901 to 2099
    daynum = (
        367 * year
        - 7 * (year + (month + 9) // 12) // 4
        + 275 * month // 9
        + day
        - 730531.5
        + greenwichtime / 24
    )

    # Mean longitude of the sun
    mean_long = daynum * 0.01720279239 + 4.894967873

    # Mean anomaly of the Sun
    mean_anom = daynum * 0.01720197034 + 6.240040768

    # Ecliptic longitude of the sun
    eclip_long = (
        mean_long
        + 0.03342305518 * sin(mean_anom)
        + 0.0003490658504 * sin(2 * mean_anom)
    )

    # Obliquity of the ecliptic
    obliquity = 0.4090877234 - 0.000000006981317008 * daynum

    # Right ascension of the sun
    rasc = atan2(cos(obliquity) * sin(eclip_long), cos(eclip_long))

    # Declination of the sun
    decl = asin(sin(obliquity) * sin(eclip_long))

    # Local sidereal time
    sidereal = 4.894961213 + 6.300388099 * daynum + rlon

    # Hour angle of the sun
    hour_ang = sidereal - rasc

    # Local elevation of the sun
    elevation = asin(sin(decl) * sin(rlat) + cos(decl) * cos(rlat) * cos(hour_ang))

    # Local azimuth of the sun
    azimuth = atan2(
        -cos(decl) * cos(rlat) * sin(hour_ang),
        sin(decl) - sin(rlat) * sin(elevation),
    )

    # Convert azimuth and elevation to degrees
    azimuth = into_range(deg(azimuth), 0, 360)
    elevation = into_range(deg(elevation), -180, 180)

    # Refraction correction (optional)
    if refraction:
        targ = rad((elevation + (10.3 / (elevation + 5.11))))
        elevation += (1.02 / tan(targ)) / 60

    # Return azimuth and elevation in degrees
    #return (round(azimuth, 2), round(elevation, 2))
    return azimuth, elevation

def solar_zenith(when, location, refraction, tilt):
    #calculate azumith (direction of sun), elevation (upward direction of sun) in degrees
    azumith, elevation = sun_position(when, location, refraction)
    zenith = 90 - elevation
    
    #conversion to radians
    azumith = math.radians(azumith)
    elevation = math.radians(elevation)
    zenith = math.radians(zenith)

    return zenith

def get_SZA(lat, lon, time_start, time_diff):
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    
    return get_altitude(lat, lon, time_obj.to_pydatetime())

def get_SZA_arr(arr):
    return get_SZA(arr[0], arr[1], str(arr[2]), float(arr[3]))


def get_SZA_array(limit, gsize, time_start, time_diff):
    # limit = [min lat, max lat, min lon, max lon]
    max_lat = int(1+round(((limit[1]-limit[0])/gsize)))
    max_lon = int(1+round(((limit[3]-limit[2])/gsize)))
    #solar_zenith = np.zeros((max_lat, max_lon))
    
    #loop through array
    """
    for i in range(0, max_lat):
        for j in range(0, max_lon):
            lat = limit[0] +  i * gsize
            lon = limit[2] + j * gsize
            
            #sza = get_SZA(lat, lon, time_start, time_diff)
            solar_zenith[i][j] = lon#sza
    """
    #apply
    #getSZA_vec = np.vectorize(get_SZA_arr)
    #t = [ [[(limit[0] +  lat * gsize), (limit[2] + lon * gsize), time_start, time_diff] for lon in range(0, max_lon)] for lat in range(0, max_lat)]
    #t = np.array(t)
    
    #solar_zenith = getSZA_vec(t)
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    solar_zenith_arr = [[solar_zenith(time_obj ,(limit[0] +  lat * gsize, limit[2] + lon * gsize), False, 0) for lon in range(0, max_lon)] for lat in range(0, max_lat)]
    #[ [get_SZA_arr([(limit[0] +  lat * gsize), (limit[2] + lon * gsize), time_start, time_diff]) for lon in range(0, max_lon)] for lat in range(0, max_lat)]
    
    return solar_zenith_arr


if __name__ == '__main__':
    limit = [-90, 90, -180, 180]
    limit = [-1, 1, -2, 2]
    gsize = 0.25
    time_start = "2020-01-01 00:00:00"
    time_diff = 30
    
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    print(time_obj)

    sza = get_SZA_array(limit, gsize, time_start, time_diff)
    print(sza)
    
