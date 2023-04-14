__author__ = "Sally Zhao"
__copyright__ = "Copyright 2023, Pyroscope"
__credits__ = ["Neil Gutkin", "Jennifer Wei", "Pawan Gupta", "Robert Levy", "Xiaohua Pan", "Zhaohui Zhang"]
__version__ = "1.0.0"
__maintainer__ = "Sally Zhao"
__email__ = "zhaosally0@gmail.com"
__status__ = "Production"
# Naming Conventions
#
# Generates metadata long_name and names depending on variable input
# Assists in statistic calculations for LEOGEO aggregation
#

# imports
import netCDF4
import numpy as np

# for variable name
sensor_references = { "MYD":"MODIS_A", "MOD":"MODIS_T", "ABI_G16":"ABI_G16", 
                     "ABI_G17":"ABI_G17", "VIIRS_SNPP":"VIIRS_SNPP", "AHI_H08":"AHI_H08",
                     "LEOGEO":"LEOGEO"}

geophys_references = { "Sensor_Zenith", "Scattering_Angle", "AOD_AllQA_550", "AOD_FilteredQA_550"}

# for long_name
sensor_references_long = { "MYD":"MODIS Aqua", "MOD":"MODIS Terra", "ABI_G16":"ABI GOES-16 (GOES-R or GOES-East)", 
                     "ABI_G17":"ABI GOES-17 (GOES-S or GOES-West)", "VIIRS_SNPP":"SNPP VIIRS",
                     "AHI_H08":"AHI Himawari-8",
                     "LEOGEO":"LEO GEO"}

geophys_references_long = { "Sensor_Zenith":"Sensor Viewing Angle", 
                      "Scattering_Angle":"", 
                      "AOD_AllQA_550":"", 
                      "AOD_FilteredQA_550":""}
statistics_references_long = {"Mean":"Mean", 
                              "STD":"Standard Deviation of",
                              "Pixels":"Number of Pixel used in calculating",
                              "TotalPixels":"Total number of level 2 pixels from all the sensors used in calculating mean gridded",
                              "SensorWeighting":"Weighting of each sensor used in calculating mean gridded"
                              }
statistics_references_long_LEOGEO = {"Mean":"Mean of gridded LEO and GEO sensors", 
                              "STD":"Standard Deviation of gridded LEO and GEO sensors ",
                              "Pixels":"Number of Pixel used in calculating",
                              "TotalPixels":"Total number of level 2 pixels from all the sensors used in calculating mean gridded LEO_GEO",
                              "SensorWeighting":"Weighting of each sensor used in calculating mean gridded",
                              "NumberOfSensors": "Number of Sensors used in calculating mean gridded LEO_GEO",
                              "SensorWeighting": "Weighting of each sensor used in calculating mean gridded LEO_GEO"
                              }

# given: geophysical variable, sensor name, statistic
# output string of nc var name
def nc_var_name(geophys_name, sensor_name, statistic = None):
    
    if "Image_Optical_Depth_Land_And_Ocean" in geophys_name:
        name = "AOD_AllQA_550"
    elif "Optical_Depth_Land_And_Ocean" in geophys_name:
        name = "AOD_FilteredQA_550"
    else:
        name = geophys_name #geophys vairable
    
    for sensor_key in sensor_references: #add appropriate sensor name
        if sensor_key in sensor_name:
            name = name + "_" + sensor_references[sensor_key] 
    
    if statistic:
        name = name + "_" + statistic
    
    return name

def nc_long_name(geophys_name, sensor_name, statistic = None, aod_long = None, filtered=False):
    name = ""
    if statistic == None:
        for sensor_key in sensor_references_long: #add appropriate sensor name
            if sensor_key in sensor_name:
                name = name + " " + sensor_references_long[sensor_key] 
        if geophys_name == "Sensor_Zenith":
            name = name + " Sensor Viewing Angle, mean for the grid"
        elif geophys_name == "Scattering_Angle":
            name = name + " Scattering Angle, mean for the grid"
        return name
    
    #AOD long name
    if sensor_name != "LEOGEO":
        for stat in statistics_references_long:
            if statistic == stat:
                name = statistics_references_long[stat]
        for sensor_key in sensor_references_long: #add appropriate sensor name
            if sensor_key in sensor_name:
                name = name + " " + sensor_references_long[sensor_key] 
                if aod_long != None:
                    name = name + " " + aod_long#" AOD at 0.55 micron (Optical_Depth_Land_And_Ocean) for both ocean (Average) (Quality flag = 1, 2, 3)  and land (corrected) (Quality flag = 3) for the grid "
    #LEOGEO
    else:
        curr_stat = ""
        for stat in statistics_references_long_LEOGEO:
            if statistic == stat:
                name = statistics_references_long_LEOGEO[stat]
                curr_stat = statistic
                if statistic =="Mean" or statistic =="STD" or statistic =="Pixels" or statistic == "TotalPixels":
                    if filtered:
                        name = name + " " + "AOD at 0.55 micron (Optical_Depth_Land_And_Ocean) for both ocean (Average) and land (corrected) with all quality data (Quality flag = 0, 1, 2, 3)"
                    else:
                        name = name + " " + "AOD at 0.55 micron (Image_Optical_Depth_Land_And_Ocean) for both ocean (Average) and land (corrected) with all quality data (Quality flag = 0, 1, 2, 3)"
        for sensor_key in sensor_references_long: #add appropriate sensor name
            if sensor_key in sensor_name:
                #name = name + " " + sensor_references_long[sensor_key] 
                if aod_long != None:
                    name = name + " " + aod_long
                if curr_stat == "SensorWeighting":
                    name = name + ", the default value is 1.0 if sensor is available and 0.0 if sensor is not available, the sensor dimension order is MODIS-T, MODIS-A, VIIRS-SNP, ABI-G16, ABI-G17, AHI-H08"
    return name

# Returns Sensor by naming convention
def get_sensor(s_name):
    for sensor_key in sensor_references: #add appropriate sensor name
        if sensor_key in s_name:
            return sensor_references[sensor_key] 
    
    return "Error"

# statistic can be Mean, STD, TotalPixels
# data = array of satellite data [[sat1 data], [sat2 data], etc]
# order is always mean, std, count
def calculate_statistic(statistic, data, data2 = None, data3 = None, scale_factor = 1):
    if statistic == "Mean":
        """
        data = np.array(data)
        num_sensors = len(data)
        lat_width = len(data[0])
        lon_width = len(data[0][0])
        data = data.flatten()
        
        data[data<=-9999.] = np.nan
        data = data.reshape(num_sensors, lat_width, lon_width)
        """
        
        data = np.array(data)
        data[data<=-9999.] = np.nan # fill_value
        
        avgtau = np.nanmean(data, axis=0 ) # [ [avgtau for modis a], [avgtau .. ]  ]
        avgtau = np.nan_to_num(avgtau, nan=-9999) # fill_value
        #print(avgtau.flatten())
        #print("ORIGINAL DATA2:", data)
        
        #avgtau = np.nanmean(np.array(data), axis=0 )
        avgtau = avgtau/scale_factor
        
        avgtau[avgtau > 5001] = -9999
        avgtau[avgtau < -101] = -9999
        
        print("AVGTAU LEOGEO TYPE: ", avgtau.dtype)
    
        return avgtau
    
    # https://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation
    if statistic == "STD":
        
        # [[std g16], [ g17], []]
        """
        mean = np.array(data)
        count = np.array(data3)
        variances = np.array(data2)**2
        
        y_bar = calculate_statistic("Mean", mean, None, count, scale_factor)
        total_count = np.sum(count, axis=0 )
        
        numerator_sum = (count - 1) * variances + count * (mean - y_bar)**2
        numerator = np.sum(numerator_sum, axis = 0)
        
        denominator = total_count - 1
        denominator[denominator == 0] = None
        
        return np.sqrt(numerator / denominator)
        """
        data = np.array(data)
        data[data<=-9999.] = np.nan
        
        stdtau = np.nanstd(data, axis = 0)
        
        return stdtau
        
    
    if statistic == "TotalPixels":
        count = np.sum(np.array(data3), axis=0 )
        return count
    
    return
    

if __name__ == '__main__':
    #testing
    data = [[[3., 2], 
             [0, -10000000]], 
            [[1, 4], 
             [1, 1]]]
    data2 = [[[10, 10], [10, 10]], [[0, 10], [10, 10]], [[10, 10], [10, 10]]]
    data3 = [[[3, 2], [0, 0]], [[1, 4], [1, -5]], [[1, 0], [0, 0]]]
    
    #print(nc_long_name("Optical_Depth_Land_And_Ocean", "LEOGEO", "TotalPixels"))
    #data = np.array([[[0,1],[2,3]], [[1,2], [3, 4]]])
    #print(np.nanmean(np.array(data), axis=0 ))
    
    #path = "/mnt/c/Users/bobgr/Desktop/NASA Spring 2023/Gridtools Package (Code, README, inputs, outputs, examples, verification)/SampleOutputs 0000-0059 01-01-2020/"
    #filename = "XAERDT_L3_MEASURES_QD_HH.20200101.0000.V0.20230130.nc"
    #L2FID = netCDF4.Dataset(path+filename,'r',format='NETCDF4')
    #L2FID.close()
    
    print(data)
    print("MEAN:", calculate_statistic("Mean", data))
    print("STD:", calculate_statistic("STD", data))