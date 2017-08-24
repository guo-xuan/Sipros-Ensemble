#!/usr/bin/python

"""
sipros_prepare_protein_database.py

sipros_prepare_protein_database.py is used created a protein database with decoy sequence
appended

Created by Xuan Guo on 02/20/2017.
Copyright (c) 2017 Xuan Guo (ORNL). Allrights reserved.
"""

## Import Python package modules
import sys, getopt, os, random

## Import Sipros package modules
import sipros_post_module
import parseconfig

## Check file exist
check_file_exist = sipros_post_module.check_file_exist

def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hi:o:c:")

    output_filename = ""
    input_filename  = ""
    config_filename = ''

    # Basic options
    for option, value in opts:
        if option in ("-h"):
            print('reverseseq.py -i input-file -o output-file -c config-file')
            sys.exit(1)
        if option in ("-i"):
            input_filename = value
        elif option in ("-o"):
            output_filename = value
        elif option in ('-c'):
            config_filename = value
        else:
            print('reverseseq.py -i input-file -o output-file -c config-file')
            sys.exit(1) 

    if input_filename == '' or output_filename == '' or config_filename == '' :
        print('reverseseq.py -i input-file -o output-file -c config-file')
        sys.exit(1)
    '''
    if (output_filename == "") :
        (inputFileNameRoot, inputFileNameExt) = os.path.splitext(input_filename)
        output_filename = inputFileNameRoot + "_CFR" + inputFileNameExt
    '''
        
    return input_filename, output_filename, config_filename


def reverse_protein_database(input_file_str, output_file_str, all_config_dict) :
    
    probability_1 = 0.5
    probability_2 = 0.5
    
    training_prefix_str = '>Rev1_'
    testing_prefix_str = '>TestRev_'
    reserved_prefix_str = '>Rev2_'
    
    if '' in all_config_dict:
        pass
    
    
    outputFile = open(output_file_str, "w")
    id_str = ""
    seq_str = ""
    seq_new_str = ""
    inputFile = open(input_file_str, "r")
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

## Parse config file
def parse_config(config_filename):

    # Save all config values to dictionary
    all_config_dict = {}    # initialize dictionay
    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)
    
    return all_config_dict

def main(argv=None):

    if argv is None:
        argv = sys.argv
        
    # parse options
    input_file_str, output_file_str, config_file_str = parse_options(argv)
    
    # parse config file
    all_config_dict = parse_config(config_file_str)
    
    # reverse sequence and save file
    reverse_protein_database(input_file_str, output_file_str, all_config_dict)
    
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




