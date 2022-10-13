# importing the required modules
import os
import argparse

#import other module files
import sat_data_input
import filter_data
import gridding
import grid_ncf
import time_conv

# error messages
INVALID_FILETYPE_MSG = "Error: Invalid file format. %s must be a .nc or .hdf file."
INVALID_PATH_MSG = "Error: Invalid file path/name. Path %s does not exist."


def grid(args):
    print("hello")
    
	
	
def main():
	# create parser object
	parser = argparse.ArgumentParser(description = "Gridding tools for satellite data transformation.") 
    command_group = parser.add_mutually_exclusive_group()

	# defining arguments for parser object
	parser.add_argument("-g", "--grid", type = str, nargs = 2,
						help = "Opens, reads, filters, and grids specified file")

	# parse the arguments from standard input
	args = parser.parse_args()
	
    print(args)

	# calling functions depending on type of argument
	if args.grid:
		grid(args)


if __name__ == "__main__":
	# calling the main function
	main()
