# due to the introduction of the time variable, we need to know the time of the sensor data
# this information is not given in the file itself, but rather the file name
# e.g. "AERDT_L2_ABI_G16.A2021008.0850.001.nc" --> 10/08/2021 8:50 AM
import datetime

# given 'yyyy/mm/dd/hr/mm' string
# converts to datetime object
def to_datetime(date_string):
    format = '%Y/%m/%d/%H/%M'
    d = datetime.datetime.strptime(date_string, format)

    return d

# function that pulls time from filename (str)
# returns datetime object
def filename_time(filename):
    time = filename.split(".")[1:3]
    date = datetime.datetime.strptime(time[0][3:], '%y%j')
    
    day_time = time[1]
    date = date.replace(hour=int(day_time[:2]), minute=int(day_time[2:]))
    
    return date

# given list of files, we need to sort them based upon their dates
# input 1d array
# output 2d array where each element represents the list of dates that fall with the inteveral

# start represents the start time (datetime object)
# end represents the end time (datetime object)
# start and end are inclusive
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
            #print(f.split("/")[-1].split(".")[0])
            if f.split("/")[-1].split(".")[0] not in d: # add to dict if not already present
                d[f.split("/")[-1].split(".")[0]] = [f]
            else:
                d[f.split("/")[-1].split(".")[0]].append(f)

        #split_files.append(list(d.values()))
        split_files.append(d)

    # missing files check
    if len(split_files) < 6:
        print("Missing sensors")
    elif len(split_files)>6:
        print("Too many sensors")
    return split_files