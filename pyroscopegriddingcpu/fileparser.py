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

import yaml
import os
from os import walk
from glob import glob


# must work for both Windows and Linux systems 
# Windows: "C:\\..."
# Linx: "/mnt/..."
def convert_path_to_linux(path):
    new_path = path.replace("\\", "/")
    new_path = new_path.replace("C:", "/mnt/c")
    
    return new_path

def convert_path_to_windows(path):
    new_path = path.replace("/mnt/c", "C:")
    new_path = new_path.replace("/", "\\")
    
    return new_path

def convert_to_other_system(path):
    if "/" in path: #linux path
        return convert_path_to_windows(path)
    else:
        return convert_path_to_linux(path)

# given path to text file contaings paths to files to read read out as list of files within
def read_file_sat_data(fileloc):
    lines = []

    with open(fileloc) as f:
        lines = [convert_path_to_linux(line.rstrip('\n')) for line in f]

    #print(lines)
    return lines

# given folder path with files inside, read out list of files within
def read_folder_sat_data(folderloc):
    #convert all to linux
    folderloc = convert_path_to_linux(folderloc)
    
    files = os.listdir(folderloc)
    #reads all files in folder and converts to path of system type
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
    print('filep======================================================')
    print(variables)
    #print(file_io)

    return grid_settings, variables, file_io

def read_directory_sat_data_old(path):
    f = []
    filelist = []
    for (dirpath, dirnames, filenames) in walk(path):
        #assume all in linux
        filenames_path = [convert_path_to_linux(dirpath)+"/"+f if convert_path_to_linux(dirpath)[-1]!="/" else convert_path_to_linux(dirpath)+f for f in filenames]
        
        #convert to system path setting
        #filenames_path = [os.path.normpath(f) for f in filenames_path]
        filelist = filelist + filenames_path
        
    return filelist

def read_directory_sat_data(path):
    suffixes = ["nc","hdf"]
    filelist=[]
    for i in suffixes:
      nclist = os.path.join(path,'**','*.'+i)
      filelist.extend( glob(nclist,recursive=True))
    print(filelist)
    return filelist
if __name__ == '__main__':
    path = "/mnt/c/Users/bobgr/Desktop/NASA Spring 2023/Gridtools Package (Code, README, inputs, outputs, examples, verification)/"
    folder = "LAADSemulation/"
    f = read_directory_sat_data(path+folder)
    f = read_folder_sat_data(path+folder)
    print(f)
