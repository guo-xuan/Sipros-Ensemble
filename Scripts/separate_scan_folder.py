# -*- coding: utf-8 -*-


import sys, getopt, warnings, os, re
from datetime import datetime, date, time


def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:n:o:",
                                    ["help",
                                     "working-dir",
                                     "file-number",
                                     "output-path"])


    # Default working dir and config file
    working_dir = "./"
    config_file = "SiprosConfig.cfg"
    filenumber  = ""
    outputpath  =""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print("-w workingdirectory/ -n filenumber -o output-path")
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-n", "--file-number"):   
            filenumber = value
        if option in ("-o", "--working-dir"):    
            outputpath = value

    if (filenumber == "") or (outputpath == "") :
        print("please specify -n and -o\n")
        sys.exit(1)                
    FT2_filename_list = get_file_list_with_ext(working_dir, ".FT2")
    MS2_filename_list = get_file_list_with_ext(working_dir, ".MS2")
    Scans_filename_list = FT2_filename_list + MS2_filename_list 

    return [Scans_filename_list, filenumber, outputpath]

## Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list

## +------+
## | Main |
## +------+
def main(argv=None):

    # try to get arguments and error handling
        if argv is None:
		argv = sys.argv
       		 # parse options
		[Scans_filename_list, filenumber, outputpath] = parse_options(argv)
        for eachFileName in Scans_filename_list :
                os.system("python separate_scans.py  "+eachFileName+" "+outputpath+" "+filenumber)



## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())




