
from pysolar.solar import get_altitude

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

def get_SZA(lat, lon, time_start, time_diff):
    time_obj = pd.to_datetime(time_start) + pd.to_timedelta(time_diff, unit="minute")
    time_obj = time_obj.replace(tzinfo=datetime.timezone.utc)
    
    return get_altitude(lat, lon, time_obj.to_pydatetime())

def get_SZA_array(limit, gsize, time_start, time_diff):
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

if __name__ == '__main__':
    limit = [-90, 90, -180, 180]
    gsize = 0.25
    time_start = "2020-01-01 00:00:00"
    time_diff = 30
    
    
    sza = get_SZA_array(limit, gsize, time_start, time_diff)
    print(sza)