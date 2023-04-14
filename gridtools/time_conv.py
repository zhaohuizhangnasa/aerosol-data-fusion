__author__ = "Sally Zhao"
__copyright__ = "Copyright 2023, Pyroscope"
__credits__ = ["Neil Gutkin", "Jennifer Wei", "Pawan Gupta", "Robert Levy", "Xiaohua Pan", "Zhaohui Zhang"]
__version__ = "1.0.0"
__maintainer__ = "Sally Zhao"
__email__ = "zhaosally0@gmail.com"
__status__ = "Production"
# Time Conversion
#
# Buckets file names according to user specifications
#

# due to the introduction of the time variable, we need to know the time of the sensor data
# this information is not given in the file itself, but rather the file name
# e.g. "AERDT_L2_ABI_G16.A2021008.0850.001.nc" --> 10/08/2021 8:50 AM
import datetime

full_satellite_list = ['ABI_G16', 'ABI_G17', 'AHI_H08', 'VIIRS_SNPP', 'MOD04', 'MYD04']

# given 'yyyy/mm/dd/hr/mm' string
# converts to datetime object
def to_datetime(date_string):
    format = '%Y/%m/%d/%H/%M'
    d = datetime.datetime.strptime(date_string, format)

    return d

# get system date
def sys_time():
    now = str(datetime.datetime.now())
    now_str = now[:10]
    now_str = now_str.replace("-", "")
    
    return str(now_str)

# function that pulls time from filename (str)
# returns datetime object
def filename_time(filename):
    #print("FILENAME: ", filename)
    time = filename.split(".")[1:3]
    date = datetime.datetime.strptime(time[0][3:], '%y%j')
    
    day_time = time[1]
    date = date.replace(hour=int(day_time[:2]), minute=int(day_time[2:]))
    
    return date

# given list of files, we need to sort them based upon their dates
# input 1d array [[], [], []]
# output 2d array where each element represents the list of dates that fall with the inteveral

# start represents the start time (datetime object)
# end represents the end time (datetime object)
# start and end are inclusive (2/1/23 - should it be???)
# interval is in minutes
# intervals will be [start, start+interval] where first index is inclusive, second not
def split_filetimes(filelist, start, end, interval):
    split_files = []
    curr = start
    end_curr = start + datetime.timedelta(minutes = interval)
    
    while (curr <= end): # keep going while current starting time interval is valid
        temp = []
        for file in filelist:
            if filename_time(file) >= curr and filename_time(file) < end_curr and filename_time(file) <= end:
                temp.append(file)
        curr = curr + datetime.timedelta(minutes = interval)
        end_curr = end_curr + datetime.timedelta(minutes = interval)
        
        split_files.append(temp)
    
    return split_files

# given a list of splitfiles (split on time)
# e.g. : [[time 0], [time1], ... [timek]]
# Additionally splits the inside lists based on common name
# Creates dictionary where key is sensor name
def split_filenames(filelist):
    split_files = []
    
    for time in filelist: # time list
        d = {} # dictionary

        for f in time:
            #retrieve sensor name from file format
            s_name = f.split("/")[-1].split(".")[0]
            
            #check to see if sensor is from list of the six measures satellites
            for measures_satellite in full_satellite_list:
                if measures_satellite in s_name:
                    s_name = measures_satellite
                    
            if s_name not in d: # add to dict if not already present
                d[s_name] = [f]
            else:
                d[s_name].append(f)

        #split_files.append(list(d.values()))
        split_files.append(d)
        
    """
    # missing files check
    if len(split_files) < 6:
        print("Missing sensors")
    elif len(split_files)>6:
        print("Too many sensors")
    """
    return split_files

#testing
if __name__ == '__main__':
    files2019 = ["XAERDT_L2_ABI_G16.A2019208.0000.001.2022260004939.nc",
                 "XAERDT_L2_ABI_G16.A2019208.0010.001.2022260004939.nc",
                 "XAERDT_L2_ABI_G16.A2019208.0020.001.2022260004940.nc",
                 "XAERDT_L2_ABI_G16.A2019208.0030.001.2022260004941.nc",
                "XAERDT_L2_ABI_G16.A2019208.0040.001.2022260004942.nc",
                "XAERDT_L2_ABI_G16.A2019208.0050.001.2022260004942.nc",
                "XAERDT_L2_ABI_G16.A2019208.0100.001.2022260004939.nc"]
    start = to_datetime("2019/07/27/00/00")
    end = to_datetime("2019/07/27/00/40")
    
    files_split_time = split_filetimes(files2019, start, end, 30)
    #print("Split on time: ", files_split_time)
    
    files_split_sensors = split_filenames(files_split_time)
    
    #print("\n\nSplit on sensors: ", files_split_sensors)
    print(files_split_sensors)
    
    for key in files_split_sensors[0]:
        print(key)
    
