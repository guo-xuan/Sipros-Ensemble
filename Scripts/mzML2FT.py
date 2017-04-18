#!/usr/bin/python

## Import Python package modules
import sys, getopt, warnings, os, re

import pymzml

def parse_options(argv):

    
    opts, args = getopt.getopt(argv[1:], "hw:o:",
                                    ["help",
                                     "working-dir",
				                     "output-dir",])


    # Default working dir and config file
    working_dir = "./"
    output_dir = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print("-w workingdirectory -o outputdirectory")
            sys.exit(1)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-o", "--output-dir"):
            output_dir = value

    mzML_filename_list = get_file_list_with_ext(working_dir, ".mzML")
    mzML_filename_list = mzML_filename_list + get_file_list_with_ext(working_dir, ".mzml")
    mzML_filename_list = mzML_filename_list + get_file_list_with_ext(working_dir, ".MZML")
    
    if (output_dir == "") :
        output_dir = working_dir

    return [mzML_filename_list, output_dir]
    
    
    
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

       # if len(file_list) == 0:
        #    print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            # die("Program exit!")
	 #   sys.exit(0)
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list


def ParseSpectrum_FT2(current_spectrum, current_file,  iParentScanNumber) :
    #print("precursors charge", current_spectrum['precursors'][0].get('charge'))
    #print("precursors mz", current_spectrum['precursors'][0].get('mz'))
    #print("S:", current_spectrum['id'])
    iScanNumber = current_spectrum['id']
    sParentMZ   = "%.5f" % current_spectrum['precursors'][0].get('mz')
    if (current_spectrum['precursors'][0].get('charge') != None ) :
        iParentCharge = current_spectrum['precursors'][0].get('charge')
        sParentMass   = "%.5f" %  (current_spectrum['precursors'][0].get('mz')  * iParentCharge)
    else :
        iParentCharge = 0
        sParentMass   = "0"
    dRetentionTime = "%.5f" %  current_spectrum['MS:1000016']
    sScanType      = "NA"
    sScanFilter    = current_spectrum['filter string']

    current_file.write("S\t"+str(iScanNumber)+"\t"+str(iScanNumber)+"\t"+sParentMZ+"\n")
    current_file.write("Z\t"+str(iParentCharge)+"\t"+ sParentMass +"\n")
    current_file.write("I\tRetentionTime\t"+str(dRetentionTime)+"\n")
    current_file.write("I\tScanType\t"+sScanType+"\n")
    current_file.write("I\tScanFilter\t"+sScanFilter+"\n")
    current_file.write("D\tParentScanNumber\t"+str(iParentScanNumber)+"\n")

    for each_fragment_mz, each_fragment_intensity in current_spectrum.centroidedPeaks :
        if each_fragment_intensity > 0 :
            seach_fragment_mz        = "%.5f" %  each_fragment_mz
            seach_fragment_intensity = "%.2f" %  each_fragment_intensity
            current_file.write(seach_fragment_mz +"\t"+seach_fragment_intensity+"\t0\t0\t0\t0\n")


def ParseSpectrum_FT1(current_spectrum, current_file) :
    #print("MS2 level:", current_spectrum['ms level'], current_spectrum['MS:1000511'])
    #print("Scan Start Time:", current_spectrum['MS:1000016'])
    #print("ScanFilter:", current_spectrum['filter string'])
    iScanNumber = current_spectrum['id']
    dRetentionTime = "%.5f" %  current_spectrum['MS:1000016']
    sScanType      = "NA"
    sScanFilter    = current_spectrum['filter string']

    current_file.write("S\t"+str(iScanNumber)+"\t"+str(iScanNumber)+"\n")
    current_file.write("I\tRetentionTime\t"+str(dRetentionTime)+"\n")
    current_file.write("I\tScanType\t"+sScanType+"\n")
    current_file.write("I\tScanFilter\t"+sScanFilter+"\n")

    for each_fragment_mz, each_fragment_intensity in current_spectrum.centroidedPeaks :
        if each_fragment_intensity > 0 :
            seach_fragment_mz        = "%.5f" %  each_fragment_mz
            seach_fragment_intensity = "%.2f" %  each_fragment_intensity
            current_file.write(seach_fragment_mz+"\t"+seach_fragment_intensity+"\t0\t0\t0\t0\n")

    return iScanNumber
    

def ConvertmzMLFile(current_mzML_filename, output_dir) :
    mzMLFileNameBase = os.path.basename(current_mzML_filename)
    (mzMLFileNameRoot, mzMLFileNameExt) = os.path.splitext(mzMLFileNameBase)
    current_FT1_filename = output_dir + os.sep + mzMLFileNameRoot + ".FT1"
    current_FT2_filename = output_dir + os.sep + mzMLFileNameRoot + ".FT2"
    current_FT1_file = file(current_FT1_filename, "w")
    current_FT2_file = file(current_FT2_filename, "w")
    msscan_list = pymzml.run.Reader(current_mzML_filename, extraAccessions=[('MS:1000515',['value'])] )
    current_FT1_file.write("H\tExtractor\tmzML2FT\n")
    current_FT2_file.write("H\tExtractor\tmzML2FT\n")
    current_FT1_file.write("H\tm/z\tIntensity\tResolution\tBaseline\tNoise\tCharge\n")
    current_FT2_file.write("H\tm/z\tIntensity\tResolution\tBaseline\tNoise\tCharge\n")
    current_FT1_file.write("H\tInstrument Model\tNA\n")
    current_FT2_file.write("H\tInstrument Model\tNA\n")

    #bFirstSpectrum = True 
    iParentScanNumber = 0
    for each_spectrum in msscan_list :
        if each_spectrum['ms level'] == 1 :
            iParentScanNumber = ParseSpectrum_FT1(each_spectrum, current_FT1_file)
        elif each_spectrum['ms level'] == 2 :
            ParseSpectrum_FT2(each_spectrum, current_FT2_file, iParentScanNumber)
            
        else :
            if (each_spectrum['ms level'] != None) :
                print("wrong ms level: "+ each_spectrum['ms level'])
            if (each_spectrum['id'] != None) and (each_spectrum['id'] != "TIC") and (each_spectrum['id'] != "tic"):
                print("the troubling Scan Number is " + each_spectrum['id'])
            elif (each_spectrum['id'] != "TIC") and (each_spectrum['id'] != "tic") :
                print("the current precursor nubmer is " + str(iParentScanNumber))
            # exit(1)

    current_FT1_file.close()
    current_FT2_file.close()

def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
        # parse options
        [mzML_filename_list, output_dir] = parse_options(argv)
    for each_mzML_filename in mzML_filename_list :
        ConvertmzMLFile(each_mzML_filename, output_dir)
        


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    main()
    

    
    
