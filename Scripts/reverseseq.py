#!/usr/bin/python


## Import Python package modules
import sys, getopt, os, random

def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hi:o:",
                                    ["help",
                             	     "input-file",
	                			     "output-file"])

    output_filename = ""
    input_filename  = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print("-i input-file, -o output-file")
            sys.exit(1)
        if option in ("-i", "--input-file"):
            input_filename = value
        if option in ("-o", "--output-file"):
            output_filename = value

    if (input_filename == "") :
        print("Please specify -i")
        sys.exit(1)
    if (output_filename == "") :
        (inputFileNameRoot, inputFileNameExt) = os.path.splitext(input_filename)
        output_filename = inputFileNameRoot + "_CFR" + inputFileNameExt
    return (input_filename, output_filename)


def ReverseSeq(inputFileName, outputFileName) :
    outputFile = open(outputFileName, "w")
    id_str = ""
    seq_str = ""
    seq_new_str = ""
    inputFile = open(inputFileName, "r")
    line_str = ""
    for line_str in inputFile:
        if line_str[0] == '>':
            if seq_str != "":
                if id_str.startswith('>Rev_') or id_str.startswith('>rev_'):
                    id_str = line_str
                    seq_str = ""
                    continue
                seq_new_str = (seq_str[::-1])
                outputFile.write(id_str)
                outputFile.write(seq_str)
                outputFile.write('\n')
                if random.random() >= 0.5:
                    outputFile.write(">Rev1_")
                else:
                    outputFile.write(">Rev2_")
                outputFile.write(id_str[1:])
                outputFile.write(seq_new_str)
                outputFile.write("\n")
            id_str = line_str
            seq_str = ""
        else:
            seq_str += line_str.strip()
    if seq_str != "":
        if id_str.startswith('>Rev_') or id_str.startswith('>rev_'):
            id_str = line_str
            seq_str = ""
        else:
            seq_new_str = (seq_str[::-1])
            outputFile.write(id_str)
            outputFile.write(seq_str)
            outputFile.write('\n')
            if random.random() >= 0.5:
                outputFile.write(">Rev1_")
            else:
                outputFile.write(">Rev2_")
            outputFile.write(id_str[1:])
            outputFile.write(seq_new_str)
            outputFile.write("\n")
        
    inputFile.close()
    outputFile.close()


def ReverseSeq_3(inputFileName, outputFileName) :
    outputFile = open(outputFileName, "w")
    id_str = ""
    seq_str = ""
    seq_new_str = ""
    inputFile = open(inputFileName, "r")
    line_str = ""
    for line_str in inputFile:
        if line_str[0] == '>':
            if seq_str != "":
                if id_str.startswith('>Rev_'):
                    id_str = line_str
                    seq_str = ""
                    continue
                seq_new_str = (seq_str[::-1])
                outputFile.write(id_str)
                outputFile.write(seq_str)
                outputFile.write('\n')
                rand_float = random.random()
                if rand_float < 0.33333333:
                    outputFile.write(">Rev_1_")
                elif rand_float < 0.66666666:
                    outputFile.write(">Rev_2_")
                else:
                    outputFile.write(">TestRev_")
                outputFile.write(id_str[1:])
                outputFile.write(seq_new_str)
                outputFile.write("\n")
            id_str = line_str
            seq_str = ""
        else:
            seq_str += line_str.strip()
    if seq_str != "":
        if id_str.startswith('>Rev_'):
            id_str = line_str
            seq_str = ""
        else:
            seq_new_str = (seq_str[::-1])
            outputFile.write(id_str)
            outputFile.write(seq_str)
            outputFile.write('\n')
            rand_float = random.random()
            if rand_float < 0.33333333:
                outputFile.write(">Rev_1_")
            elif rand_float < 0.66666666:
                outputFile.write(">Rev_2_")
            else:
                outputFile.write(">TestRev_")
            outputFile.write(id_str[1:])
            outputFile.write(seq_new_str)
            outputFile.write("\n")
        
    inputFile.close()
    outputFile.close()


def main(argv=None):

    # try to get arguments and error handling
    if argv is None:
        argv = sys.argv
        # parse options
    (inputFileName, outputFileName) = parse_options(argv)
    # ReverseSeq_3(inputFileName, outputFileName)
    ReverseSeq(inputFileName, outputFileName)
    
    print('Done.')    
    
def main2(argv=None):
    folder_str = "/media/xgo/Seagate/Proteomics/Data/Zhou/Pete/2014_database/"
    output_str = "/media/xgo/Seagate/Proteomics/Data/Zhou/Pete/2014_database_fwr_rev/"
    for file_str in os.listdir(folder_str):
        if file_str.endswith('.faa'):
            ReverseSeq(folder_str+file_str, output_str+file_str[:-4]+'_fwr_rev.faa')
            print('.')
    
    print('Done.')


## If this program runs as standalone, then go to main.
if __name__ == "__main__":
    sys.exit(main())




