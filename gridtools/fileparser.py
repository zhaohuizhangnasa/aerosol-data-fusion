__author__ = "Sally Zhao"
__copyright__ = "Copyright 2023, Pyroscope"
__credits__ = ["Neil Gutkin", "Jennifer Wei", "Pawan Gupta", "Robert Levy", "Xiaohua Pan", "Zhaohui Zhang"]
__version__ = "1.0.0"
__maintainer__ = "Sally Zhao"
__email__ = "zhaosally0@gmail.com"
__status__ = "Production"
# FileParser
#
# Aids in parsing file names and paths
# Controls directory read in
#

import time_conv
import yaml
import os
from os import walk

# given path to text file contaings paths to files to read read out as list of files within
def read_file_sat_data(fileloc):
    lines = []

    with open(fileloc) as f:
        lines = [line.rstrip('\n') for line in f]

    #print(lines)
    return lines

# given folder path with files inside, read out list of files within
def read_folder_sat_data(folderloc):
    files = os.listdir(folderloc)
    lines = [folderloc+"/"+f for f in files]
    
    return lines

# netCDF_multi_time(filelist, gsize, geo_list, phy_list, output, time_interval, time_start, time_end)
def read_in_commands(fileloc):
    cmds = []

    with open(fileloc) as f:
        cmds = [line.rstrip('\n') for line in f]

    filelist = read_file_sat_data(cmds[1])
    gsize = float(cmds[2])
    time_interval = int(cmds[3])
    start = time_conv.to_datetime(cmds[4])
    end = time_conv.to_datetime(cmds[5])
    output = cmds[6]

    geo_list = cmds[7].split(" ")
    phy_list = cmds[8].split(" ")

    return filelist, gsize, time_interval, start, end, output, geo_list, phy_list

def read_config(yfile):
    #print("FILE: ", yfile)
    #print("OS: ", os.getcwd())
    with open(yfile, "r") as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    grid_settings = data["grid_settings"]
    variables = data["variables"]
    file_io = data["file_io"]

    #print(grid_settings)
    #print(variables)
    #print(file_io)

    return grid_settings, variables, file_io

def read_directory_sat_data(path):
    f = []
    filelist = []
    for (dirpath, dirnames, filenames) in walk(path):
        #print(dirpath)
        #print(dirnames)
        #print(filenames)
        filenames_path = [dirpath+"/"+f if dirpath[-1]!="/" else dirpath+f for f in filenames]
        filelist = filelist + filenames_path
        
    return filelist

if __name__ == '__main__':
    path = "/mnt/c/Users/bobgr/Desktop/NASA Spring 2023/Gridtools Package (Code, README, inputs, outputs, examples, verification)/"
    folder = "LAADSemulation/"
    f = read_directory_sat_data(path+folder)
    f = read_folder_sat_data(path+folder)
    print(f)
