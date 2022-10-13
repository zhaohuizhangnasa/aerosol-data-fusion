# imports
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
                              "STD":"Standard Deviation in",
                              "Pixels":"Number of Pixel used in calculating",
                              "TotalPixels":"Total number of level 2 pixels from all the sensors used in calculating mean gridded",
                              "SensorWeighting":"Weighting of each sensor used in calculating mean gridded"
                              }
statistics_references_long_LEOGEO = {"Mean":"Mean of gridded LEO and GEO sensors", 
                              "STD":"Standard Deviation of gridded LEO and GEO sensors ",
                              "Pixels":"Number of Pixel used in calculating",
                              "TotalPixels":"Total number of level 2 pixels from all the sensors used in calculating mean gridded mean gridded LEO_GEO",
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

def nc_long_name(geophys_name, sensor_name, statistic = None, aod_long = None):
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
                name = name + " " + aod_long#" AOD at 0.55 micron (Optical_Depth_Land_And_Ocean) for both ocean (Average) (Quality flag = 1, 2, 3)  and land (corrected) (Quality flag = 3) for the grid "
    #LEOGEO
    else:
        curr_stat = ""
        for stat in statistics_references_long:
            if statistic == stat:
                name = statistics_references_long_LEOGEO[stat]
                curr_stat = statistic
        for sensor_key in sensor_references_long: #add appropriate sensor name
            if sensor_key in sensor_name:
                #name = name + " " + sensor_references_long[sensor_key] 
                name = name + " " + aod_long#"
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
def calculate_statistic(statistic, data, data2 = None, data3 = None):
    if statistic == "Mean":
        """
        count = np.array(data3)
        
        numerator = np.sum(np.array(data) * count, axis=0)
        denominator = np.sum(count, axis=0)
        denominator[denominator == 0] = None
        avgtau = numerator / denominator"""
        avgtau = np.nanmean(np.array(data), axis=0 )
        return avgtau
    
    # https://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation
    if statistic == "STD":
        mean = np.array(data)
        count = np.array(data3)
        variances = np.array(data2)**2
        
        y_bar = calculate_statistic("Mean", mean, None, count)
        total_count = np.sum(count, axis=0 )
        
        numerator_sum = (count - 1) * variances + count * (mean - y_bar)**2
        numerator = np.sum(numerator_sum, axis = 0)
        
        denominator = total_count - 1
        denominator[denominator == 0] = None
        
        return np.sqrt(numerator / denominator)
    
    if statistic == "TotalPixels":
        count = np.sum(np.array(data3), axis=0 )
        return count
    
    return
    

if __name__ == '__main__':
    #testing
    data = [[[3, 2], [0, 0]], [[1, 4], [1, -5]], [[1, 0], [0, 0]]]
    data2 = [[[10, 10], [10, 10]], [[0, 10], [10, 10]], [[10, 10], [10, 10]]]
    data3 = [[[3, 2], [0, 0]], [[1, 4], [1, -5]], [[1, 0], [0, 0]]]
    
    print(nc_long_name("Optical_Depth_Land_And_Ocean", "LEOGEO", "NumberOfSensors", "asdfsdf"))