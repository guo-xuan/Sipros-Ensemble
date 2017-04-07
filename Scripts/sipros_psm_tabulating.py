'''
Created on Jun 28, 2016

@author: Xuan Guo
'''

import sys, os
import getopt
import csv
from collections import defaultdict
from datetime import datetime, date, time
import re

# # Import Sipros package modules
import sipros_post_module
import parseconfig
import sipros_post_processing
from multiprocessing import Queue, cpu_count

## Returns the current time in a nice format
curr_time = sipros_post_module.curr_time
## Format time as a pretty string
format_time = sipros_post_module.format_time

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_psm_tabulating.py [options]

Inputs:
    spe2psm directory
    sipros configuration file
    output directory

Options:
    -h/--help
    -v/--version
    -i/--input-folder ./path    
    -c/--configuration
    -o/--output-folder ./path

Outputs:
    output PSM table
'''

spe2psm_file_ext = 'Spe2Pep.txt'
# # Get file(s) list in working dir with specific file extension
get_file_list_with_ext = sipros_post_module.get_file_list_with_ext
# # Class for sipros fields object
SpectrumFields = sipros_post_module.SpectrumFields
# #
PsmFields = sipros_post_module.PsmFields
# # Class for ignoring comments '#' in sipros file
CommentedFile = sipros_post_module.CommentedFile
# # check_file_exist
check_file_exist = sipros_post_module.check_file_exist

# # Decoy Reverse Forward protein
def protein_type(protein_sequence):
    sProteins = protein_sequence
    asProteins = sProteins.split(',')
    for sProtein in asProteins:
        if not (sProtein.startswith('Rev_') or sProtein.startswith('Dec_')):
            return 1
    return -1

pep_iden_str = '[Peptide_Identification]'
search_name_str = 'Search_Name'
FASTA_Database_str = 'FASTA_Database'
Maximum_Missed_Cleavages_str = 'Maximum_Missed_Cleavages'
Cleave_After_Residues_str = 'Cleave_After_Residues'

# # Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}  # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    # Save all config values to dictionary
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)
    (modification_dict, element_modification_list_dict) = parseconfig.getModificationDictionary(all_config_dict)
    # valiables were defined in global
    # pep_iden_str      = '[Protein_Identification]'
    # cleave_after_str  = 'Decoy_Prefix'
    # cleave_before_str = 'FDR_Filtering'
    # FDR_threshold_str = 'FDR_Threshold'

    # only save protein_identification config info to config_dict
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + search_name_str):
            config_dict[search_name_str] = value
        elif key == (pep_iden_str + FASTA_Database_str):
            config_dict[FASTA_Database_str] = value
        elif key == (pep_iden_str + Maximum_Missed_Cleavages_str):
            config_dict[Maximum_Missed_Cleavages_str] = value
        elif key == (pep_iden_str + Cleave_After_Residues_str):
            config_dict[Cleave_After_Residues_str] = value
        else:
            continue

    # return config dictionary
    return (config_dict, modification_dict, element_modification_list_dict)


# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:x",
                                    ["help",
                                     "version",
                                     "input-folder",
                                     "configuration",
                                     "output-folder",
                                     "pepxml"])

    # Default working dir and config file
    input_folder = ""
    sConfig = ""
    output_folder = ""
    pepxml_output = False

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print help_message
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print "sipros_peptides_filtering.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input-folder"):
            input_folder = value
        if option in ("-c", "--configuration"):
            sConfig = value
        if option in ("-o", "--output-folder"):
            output_folder = value
        if option in ("-x", "--pepxml"):
            pepxml_output = True

    if input_folder == "" or sConfig == "" or output_folder == "":
        print help_message
        sys.exit(0)

    return (input_folder, sConfig, output_folder, pepxml_output)

# # Read spe2psm files
def read_sipros_files(spe2psm_file_list):
    # save data for searching cut-off threshold
    sipros_psm_data = defaultdict(list)  # initialize dict of list
    list_file_name = []
    psm_Filename = ""
    psm_ScanNumber = ""
    psm_MeasuredParentMass = ""
    psm_Charge = ""
    psm_Mvh = []
    psm_Wdp = []
    psm_Xcorr = []
    psm_id = ""
    global iNumScoreEqual
    iNumScoreEqual = 0
    for _file_idx, sipros_file in enumerate(spe2psm_file_list):
        # read line with csv
        sipros_reader = csv.reader(CommentedFile(open(sipros_file, 'rb')),
                                   delimiter='\t')
        # skip header
        _headline = sipros_reader.next()
        _headline = sipros_reader.next()
        for _line_idx, line in enumerate(sipros_reader):
            if line[0] == '+':
                sipros_obj = SpectrumFields._make(line)
                psm_Filename = sipros_obj.Filename.strip()
                psm_ScanNumber = sipros_obj.ScanNumber.strip()
                psm_MeasuredParentMass = sipros_obj.MeasuredParentMass
                psm_Charge = sipros_obj.ParentCharge
                psm_Mvh = []
                psm_Wdp = []
                psm_Xcorr = []
                psm_id = 0
                if psm_Filename in list_file_name:
                    psm_id = list_file_name.index(psm_Filename)
                else:
                    list_file_name.append(psm_Filename)
                    psm_id = list_file_name.index(psm_Filename)
                psm_ID = str(psm_id) + '_+_' + psm_ScanNumber.strip()
                if psm_ID not in sipros_psm_data:
                    # save psm_out_list to the sipros_psm_data list
                    psm_data_list = [psm_id,
                                 psm_ScanNumber,
                                 psm_MeasuredParentMass,
                                 psm_Charge,
                                 psm_Mvh,
                                 psm_Wdp,
                                 psm_Xcorr,
                                 0]
                    # save sipros_psm_data
                    sipros_psm_data[psm_ID] = psm_data_list
            elif line[0] == '*':
                sipros_psm_obj = PsmFields._make(line)
                psm_IdentifiedPeptide = sipros_psm_obj.IdentifiedPeptide.strip()[1:-1]
                psm_OriginalPeptide = sipros_psm_obj.OriginalPeptide.strip()[1:-1]
                psm_CalculatedParentMass = sipros_psm_obj.CalculatedParentMass.strip()
                psm_MVH = float(sipros_psm_obj.MVH.strip())
                psm_WDP = float(sipros_psm_obj.WDP.strip())
                psm_Xcorr = float(sipros_psm_obj.Xcorr.strip())
                psm_ProteinNames = sipros_psm_obj.ProteinNames.strip()[1:-1]
                psm_data_pep = [psm_IdentifiedPeptide,
                                psm_OriginalPeptide,
                                psm_CalculatedParentMass,
                                psm_MVH,
                                psm_WDP,
                                psm_Xcorr,
                                psm_ProteinNames]
                # set the best MVH
                psm_data_list = sipros_psm_data[psm_ID]
                if not psm_data_list[4]:
                    psm_data_list[4] = psm_data_pep
                else:
                    if psm_data_pep[3] > psm_data_list[4][3]:
                        psm_data_list[4] = psm_data_pep
                    elif psm_data_pep[3] == psm_data_list[4][3]:
                        if psm_data_pep[0] == psm_data_list[4][0]:
                            #pass
                            psm_data_list[4][6] += ',' + psm_data_pep[6]
                        iNumScoreEqual += 1
                # set the best WDP
                if not psm_data_list[5]:
                    psm_data_list[5] = psm_data_pep
                else:
                    if psm_data_pep[4] > psm_data_list[5][4]:
                        psm_data_list[5] = psm_data_pep
                    elif psm_data_pep[4] == psm_data_list[5][4]:
                        if psm_data_pep[0] == psm_data_list[5][0]:
                            #pass
                            psm_data_list[5][6] += ',' + psm_data_pep[6]
                        iNumScoreEqual += 1
                # set the best Xcorr
                if not psm_data_list[6]:
                    psm_data_list[6] = psm_data_pep
                else:
                    if psm_data_pep[5] > psm_data_list[6][5]:
                        psm_data_list[6] = psm_data_pep
                    elif psm_data_pep[5] == psm_data_list[6][5]:
                        if psm_data_pep[0] == psm_data_list[6][0]:
                            #pass
                            psm_data_list[6][6] += ',' + psm_data_pep[6]
                        iNumScoreEqual += 1

    return sipros_psm_data

class Scores:
    lTargetScores = []
    lDecoyScores = []

    def __init__(self):
        pass

    def addTargetScores(self, fScore):
        self.lTargetScores.append(fScore)

    def addDecoyScores(self, fScore):
        self.lDecoyScores.append(fScore)

    def writeIntoFile(self, sFileName):
        with open(sFileName, 'w') as f:
            iLen = 0
            if len(self.lTargetScores) > len(self.lDecoyScores):
                iLen = len(self.lTargetScores)
            else:
                iLen = len(self.lDecoyScores)
            for i in range(iLen):
                if i < len(self.lTargetScores):
                    f.write(str(self.lTargetScores[i]))
                f.write('\t')
                if i < len(self.lDecoyScores):
                    f.write(str(self.lDecoyScores[i]))
                f.write('\n')

def report_output(sipros_psm_data):
    iNumMvh = 0
    iTargeMvh = 0
    iDecoyMvh = 0
    iNumWdp = 0
    iTargeWdp = 0
    iDecoyWdp = 0
    iNumXcorr = 0
    iTargeXcorr = 0
    iDecoyXcorr = 0
    iNumMvhWdp = 0
    iTargeMvhWdp = 0
    iDecoyMvhWdp = 0
    iNumMvhXcorr = 0
    iTargeMvhXcorr = 0
    iDecoyMvhXcorr = 0
    iNumWdpXcorr = 0
    iTargeWdpXcorr = 0
    iDecoyWdpXcorr = 0
    iNumMvhWdpXcorr = 0
    iTargeMvhWdpXcorr = 0
    iDecoyMvhWdpXcorr = 0
    lMvhScores = []
    lWdpScores = []
    lXcorrScores = []
    iIndexMvhWdpXcorr = 3
    iIndexMvhWdp = 2
    iIndexMvhXcorr = 1
    iIndexWdpMvh = 2
    iIndexWdpXcorr = 1
    iIndexXcorrMvh = 2
    iIndexXcorrWdp = 1
    for _i in range(4):
        lMvhScores.append(Scores())
        lWdpScores.append(Scores())
        lXcorrScores.append(Scores())
    lastPsm = 0
    for _psm_ID, psm_data_list in sipros_psm_data.iteritems():
        if lastPsm < int(psm_data_list[1]):
            lastPsm = int(psm_data_list[1])
        # MVH == WDP == Xcorr
        if psm_data_list[4][1] == psm_data_list[5][1] and psm_data_list[4][1] == psm_data_list[6][1]:
            iNumMvhWdpXcorr += 1
            psm_data_list[7] = 4
            if protein_type(psm_data_list[4][6]) == 1:
                iTargeMvhWdpXcorr += 1
                lMvhScores[iIndexMvhWdpXcorr].addTargetScores(psm_data_list[4][3])
                lWdpScores[iIndexMvhWdpXcorr].addTargetScores(psm_data_list[5][4])
                lXcorrScores[iIndexMvhWdpXcorr].addTargetScores(psm_data_list[6][5])
            else:
                iDecoyMvhWdpXcorr += 1
                lMvhScores[iIndexMvhWdpXcorr].addDecoyScores(psm_data_list[4][3])
                lWdpScores[iIndexMvhWdpXcorr].addDecoyScores(psm_data_list[5][4])
                lXcorrScores[iIndexMvhWdpXcorr].addDecoyScores(psm_data_list[6][5])
        # MVH == WDP != Xcorr
        elif psm_data_list[4][1] == psm_data_list[5][1]:
            iNumMvhWdp += 1
            iNumXcorr += 1
            psm_data_list[7] = 4
            if protein_type(psm_data_list[4][6]) == 1:
                iTargeMvhWdp += 1
                lMvhScores[iIndexMvhWdp].addTargetScores(psm_data_list[4][3])
                lWdpScores[iIndexWdpMvh].addTargetScores(psm_data_list[5][4])
            else:
                iDecoyMvhWdp += 1
                lMvhScores[iIndexMvhWdp].addDecoyScores(psm_data_list[4][3])
                lWdpScores[iIndexWdpMvh].addDecoyScores(psm_data_list[5][4])
            if protein_type(psm_data_list[6][6]) == 1:
                iTargeXcorr += 1
                lXcorrScores[0].addTargetScores(psm_data_list[6][5])
            else:
                iDecoyXcorr += 1
                lXcorrScores[0].addDecoyScores(psm_data_list[6][5])
        # MVH == Xcorr != WDP
        elif psm_data_list[4][1] == psm_data_list[6][1]:
            iNumMvhXcorr += 1
            iNumWdp += 1
            psm_data_list[7] = 4
            if protein_type(psm_data_list[4][6]) == 1:
                iTargeMvhXcorr += 1
                lMvhScores[iIndexMvhXcorr].addTargetScores(psm_data_list[4][3])
                lXcorrScores[iIndexXcorrMvh].addTargetScores(psm_data_list[6][5])
            else:
                iDecoyMvhXcorr += 1
                lMvhScores[iIndexMvhXcorr].addDecoyScores(psm_data_list[4][3])
                lXcorrScores[iIndexXcorrMvh].addDecoyScores(psm_data_list[6][5])
            if protein_type(psm_data_list[5][6]) == 1:
                iTargeWdp += 1
                lWdpScores[0].addTargetScores(psm_data_list[5][4])
            else:
                iDecoyWdp += 1
                lWdpScores[0].addDecoyScores(psm_data_list[5][4])
        # MVH != WDP == Xcorr
        elif psm_data_list[5][1] == psm_data_list[6][1]:
            iNumWdpXcorr += 1
            iNumMvh += 1
            psm_data_list[7] = 5
            if protein_type(psm_data_list[4][6]) == 1:
                iTargeMvh += 1
                lWdpScores[iIndexWdpXcorr].addTargetScores(psm_data_list[5][4])
                lXcorrScores[iIndexXcorrWdp].addTargetScores(psm_data_list[6][5])
            else:
                iDecoyMvh += 1
                lWdpScores[iIndexWdpXcorr].addDecoyScores(psm_data_list[5][4])
                lXcorrScores[iIndexXcorrWdp].addDecoyScores(psm_data_list[6][5])
            if protein_type(psm_data_list[5][6]) == 1:
                iTargeWdpXcorr += 1
                lMvhScores[0].addTargetScores(psm_data_list[4][3])
            else:
                iDecoyWdpXcorr += 1
                lMvhScores[0].addDecoyScores(psm_data_list[4][3])
        else:
            iNumMvh += 1
            iNumWdp += 1
            iNumXcorr += 1
            if protein_type(psm_data_list[4][6]) == 1:
                iTargeMvh += 1
                lMvhScores[0].addTargetScores(psm_data_list[4][3])
            else:
                iDecoyMvh += 1
                lMvhScores[0].addDecoyScores(psm_data_list[4][3])
            if protein_type(psm_data_list[5][6]) == 1:
                iTargeWdp += 1
                lWdpScores[0].addTargetScores(psm_data_list[5][4])
            else:
                iDecoyWdp += 1
                lWdpScores[0].addDecoyScores(psm_data_list[5][4])
            if protein_type(psm_data_list[6][6]) == 1:
                iTargeXcorr += 1
                lXcorrScores[0].addTargetScores(psm_data_list[6][5])
            else:
                iDecoyXcorr += 1
                lXcorrScores[0].addDecoyScores(psm_data_list[6][5])
    print "Mvh:\t" + str(iNumMvh) + "\tTarget/Decoy\t" + str(iTargeMvh) + "\t" + str(iDecoyMvh)
    print "Xcorr:\t" + str(iNumXcorr) + "\tTarget/Decoy\t" + str(iTargeXcorr) + "\t" + str(iDecoyXcorr)
    print "Wdp:\t" + str(iNumWdp) + "\tTarget/Decoy\t" + str(iTargeWdp) + "\t" + str(iDecoyWdp)
    print "Mvh+Xcorr:\t" + str(iNumMvhXcorr) + "\tTarget/Decoy\t" + str(iTargeMvhXcorr) + "\t" + str(iDecoyMvhXcorr)
    print "Mvh+Wdp:\t" + str(iNumMvhWdp) + "\tTarget/Decoy\t" + str(iTargeMvhWdp) + "\t" + str(iDecoyMvhWdp)
    print "Wdp+Xcorr:\t" + str(iNumWdpXcorr) + "\tTarget/Decoy\t" + str(iTargeWdpXcorr) + "\t" + str(iDecoyWdpXcorr)
    print "Mvh+Wdp+Xcorr:\t" + str(iNumMvhWdpXcorr) + "\tTarget/Decoy\t" + str(iTargeMvhWdpXcorr) + "\t" + str(iDecoyMvhWdpXcorr)
    lMvhScores[0].writeIntoFile('MvhOnly')
    lMvhScores[iIndexMvhWdp].writeIntoFile('MvhWdp')
    lMvhScores[iIndexMvhXcorr].writeIntoFile('MvhXcorr')
    lMvhScores[iIndexMvhWdpXcorr].writeIntoFile('MvhWdpXcorr')
    lWdpScores[0].writeIntoFile('WdpOnly')
    lWdpScores[iIndexWdpMvh].writeIntoFile('WdpMvh')
    lWdpScores[iIndexWdpXcorr].writeIntoFile('WdpXcorr')
    lWdpScores[iIndexMvhWdpXcorr].writeIntoFile('WdpMvhXcorr')
    lXcorrScores[0].writeIntoFile('MvhOnly')
    lXcorrScores[iIndexXcorrWdp].writeIntoFile('XcorrWdp')
    lXcorrScores[iIndexXcorrMvh].writeIntoFile('XcorrMvh')
    lXcorrScores[iIndexMvhWdpXcorr].writeIntoFile('XcorrMvhWdp')
    
    print '\nid:\t%s' % lastPsm
    return sipros_psm_data

def refine_proteins(file_str):
    pep_pro_dict = {}
    psm_list = []
    with open(file_str, 'r') as fr:
        header = fr.next()
        for line_str in fr:
            words = line_str.strip().split()
            psm_list.append(words)
            pep_str = words[7]
            pro_str = words[12]
            if pep_str in pep_pro_dict:
                pro_list = pep_pro_dict[pep_str]
                pro_words = pro_str[1:-1].split(',')
                for sProtein in pro_words:
                    if not sProtein in pro_list:
                        pro_list.append(sProtein)
            else:
                pro_list = []
                pro_words = pro_str[1:-1].split(',')
                pro_list.extend(pro_words)
                pep_pro_dict[pep_str] = pro_list
    
    with open(file_str, 'w') as fw:
        fw.write(header)
        for psm in psm_list:
            if psm[7] in pep_pro_dict:
                pro_list = pep_pro_dict[psm[7]]
                pro_str = '{' + ','.join(pro_list) + '}'
                psm[12] = pro_str
                fw.write('\t'.join(psm))
                fw.write('\n')
            else:
                print 'error'

# # simple count the agreements
def main2(argv=None):

    if argv is None:
        argv = sys.argv
    # parse options
    (input_folder, _sConfig, _output_folder) = parse_options(argv)
    # Get sipros output file(s) in working directory
    spe2psm_file_list = get_file_list_with_ext(input_folder, spe2psm_file_ext)
    # Read sipros output files and get data
    sys.stderr.write('[Step 2] Load sipros output file(s):        Running -> ')
    sipros_psm_data = read_sipros_files(spe2psm_file_list)
    sys.stderr.write('Done!\n')
    sys.stderr.write("Number Equal Scores:\t" + str(iNumScoreEqual) + "\n")
    # Report output
    sys.stderr.write('[Step 4] Report .psm.txt and .pep.txt:      Running -> \n')
    # Report output files
    report_output(sipros_psm_data)
    sys.stderr.write('Done!\n')

Spe2PepReader = sipros_post_module.Spe2PepReader
RankPsm = sipros_post_module.RankPsm
writePsm = sipros_post_module.writePsm
#writePin = sipros_post_module.writePin
#write_PepXML = sipros_post_module.write_PepXML

## Get base_out filename
get_base_out = sipros_post_module.get_base_out
PsmFields4 = sipros_post_processing.PsmFields4

class scan_xml:
    def __init__(self, line_str):
        PsmFields_obj = PsmFields4._make(line_str)
        self.filename = PsmFields_obj.FileName
        self.scannumber = PsmFields_obj.ScanNumber
        self.charge = PsmFields_obj.ParentCharge
        self.measuredmass = PsmFields_obj.MeasuredParentMass
        self.peplist = []
        self.mvhlist = []
        self.xcorrlist = []
        self.wdplist = []
        pep = pep_xml(PsmFields_obj)
        self.peplist.append(pep)
        self.list_list = []
        
    def add(self, line_str):
        PsmFields_obj = PsmFields4._make(line_str)
        pep = pep_xml(PsmFields_obj)
        self.peplist.append(pep)
        
    def score_process(self):
        self.mvhlist = sorted(self.peplist, key=lambda pep: (pep.scorelist[0]), reverse=True)
        self.xcorrlist = sorted(self.peplist, key=lambda pep: (pep.scorelist[1]), reverse=True)
        self.wdplist = sorted(self.peplist, key=lambda pep: (pep.scorelist[2]), reverse=True)
        self.list_list = [self.mvhlist, self.xcorrlist, self.wdplist]
        

class pep_xml:
    def __init__(self, PsmFields_obj):
        self.originalPep = PsmFields_obj.OriginalPeptide
        self.identifiedPep = PsmFields_obj.IdentifiedPeptide
        self.proteinlist = re.sub('[{}]', '', PsmFields_obj.ProteinNames).split(',')
        self.calculatedmass = PsmFields_obj.CalculatedParentMass
        self.scorelist = [float(PsmFields_obj.MVH), float(PsmFields_obj.Xcorr), float(PsmFields_obj.WDP)]
        self.scorediff = [float(PsmFields_obj.DiffRS1), float(PsmFields_obj.DiffRS2), float(PsmFields_obj.DiffRS3)]

fNeutronMass = 1.00867108694132 # it is Neutron mass

def get_mass_diff(measured_mass, calculated_mass):
        MassDiffOriginal = measured_mass - calculated_mass
        MassDiff = MassDiffOriginal
        for i in range(-4, 4):
            if abs(MassDiffOriginal - i*fNeutronMass) < abs(MassDiff):
                MassDiff = MassDiffOriginal - i*fNeutronMass
        return MassDiff

def get_num_missed_cleavages(peptide_str, cleave_residues_str):
    num = 0
    for c in cleave_residues_str:
        num += peptide_str.count(c)

    return str(num - 1)

Cleave_After_Residues_str = 'Cleave_After_Residues'

def get_modification_info(peptide_str, modification_label_dict):
    modification_dict = {}
    for key, value in modification_label_dict.iteritems():
        if key.isalpha():
            continue
        beg = -1
        beg = peptide_str.find(key, beg + 1)
        while beg != -1:
            modification_dict[str(beg - len(modification_dict))] = value + modification_label_dict[peptide_str[beg-1:beg]]
            beg = peptide_str.find(key, beg + 1)
    return modification_dict

def writePepxml(filename_str, config_dict, modification_dict, element_modification_list_dict):

    # start reading files
    filename_tab_str = ""
    scannumber_tab_str = ""
    psm_obj = None
    psm_list = []
    with open(filename_str, 'r') as f:
        f.next()
        for s in f:
            l = s.split('\t')
            if filename_tab_str != l[0] or scannumber_tab_str != l[1]:
                if scannumber_tab_str != "":
                    psm_obj.score_process()
                    psm_list.append(psm_obj)
                    psm_obj = None
            if psm_obj is None:
                psm_obj = scan_xml(l)
            else:
                psm_obj.add(l)
            filename_tab_str = psm_obj.filename
            scannumber_tab_str = psm_obj.scannumber            
    if psm_obj is not None:
        psm_obj.score_process()
        psm_list.append(psm_obj)
        
    writePepxmlSingle(filename_str, config_dict, modification_dict, element_modification_list_dict, psm_list, 0)
    writePepxmlSingle(filename_str, config_dict, modification_dict, element_modification_list_dict, psm_list, 1)
    writePepxmlSingle(filename_str, config_dict, modification_dict, element_modification_list_dict, psm_list, 2)
    
    
def writePepxmlSingle(filename_str, config_dict, modification_dict, element_modification_list_dict, psm_list, op):
    score_list_str = ['mvh', 'xcorr', 'wdp']
    deltascore_str = 'deltascore'
    cleave_residues_str = config_dict[Cleave_After_Residues_str]
    
    _temp = __import__('lxml.etree', globals(), locals(), ['ElementTree'], -1)
    ElementTree = _temp.ElementTree
    _temp = __import__('lxml.etree', globals(), locals(), ['Element'], -1)
    Element = _temp.Element
    _temp = __import__('lxml.etree', globals(), locals(), ['SubElement'], -1)
    SubElement = _temp.SubElement
    
    iIndexUnique = 1
    
    filename, file_extension = os.path.splitext(filename_str)
    xmlns = "http://regis-web.systemsbiology.net/pepXML"
    xsi = "http://www.w3.org/2001/XMLSchema-instance"
    schemaLocation = "http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd"
    root = Element("{" + xmlns + "}msms_pipeline_analysis",
                   nsmap={'xsi':xsi, None:xmlns},
                   attrib={"{" + xsi + "}schemaLocation"  : schemaLocation})
    root.set('date', datetime.now().isoformat())
    root.set('summary_xml', (filename + '.pepXML'))

    analysis_summary = SubElement(root, 'analysis_summary')
    analysis_summary.set('analysis', "Sipros Ensemble "+score_list_str[op])
    analysis_summary.set('version', "V1.0")
    analysis_summary.set('time', (datetime.now().isoformat()))
    
    msms_run_summary = SubElement(root, 'msms_run_summary')
    msms_run_summary.set('base_name', filename)
    msms_run_summary.set('raw_data_type', "raw")
    msms_run_summary.set('raw_data', "ms2")

    sample_enzyme = SubElement(msms_run_summary, 'sample_enzyme')
    sample_enzyme.set('name', "Trypsin/P")
    sample_enzyme.set('independent', "false")
    sample_enzyme.set('fidelity', "specific")
    
    specificity = SubElement(sample_enzyme, 'specificity')
    specificity.set('sense', "C")
    specificity.set('cut', "KR")
    specificity.set('no_cut', "")
    specificity.set('min_spacing', "1")

    search_summary = SubElement(msms_run_summary, "search_summary")
    search_summary.set('base_name', filename)
    search_summary.set('search_engine', "Sipros Ensemble "+score_list_str[op])
    search_summary.set('precursor_mass_type', 'monoisotopic')
    search_summary.set('fragment_mass_type', 'monoisotopic')
    search_summary.set('out_data_type', '')
    search_summary.set('out_data', '')
    
    search_database = SubElement(search_summary, 'search_database')
    (path,file) = os.path.split(config_dict[FASTA_Database_str])
    search_database.set('local_path', path)
    search_database.set('database_name', file)
    search_database.set('type', 'AA')

    enzymatic_search_constraint = SubElement(search_summary, 'enzymatic_search_constraint')
    enzymatic_search_constraint.set('enzyme', 'Trypsin/P')
    enzymatic_search_constraint.set('max_num_internal_cleavages', config_dict[Maximum_Missed_Cleavages_str])
    enzymatic_search_constraint.set('min_number_termini', '2')
    
    for e_key, e_value in element_modification_list_dict.iteritems():
        for e_2_value in e_value:
            aminoacid_modification = SubElement(search_summary, 'aminoacid_modification')
            aminoacid_modification.set('aminoacid', e_key)
            aminoacid_modification.set('massdiff', str(modification_dict[e_2_value]))
            aminoacid_modification.set('mass', str(modification_dict[e_2_value] + modification_dict[e_key]))
            aminoacid_modification.set('variable', e_2_value)

    # create scan 
    for psm_obj in psm_list:
        # query results
        spectrum_query = SubElement(msms_run_summary, 'spectrum_query')
        filename_tab_str = psm_obj.filename
        ScanNumber_str = psm_obj.scannumber
        ParentCharge_str = psm_obj.charge
        spectrum_query.set('spectrum', filename_tab_str + '.' + ScanNumber_str + '.' + ScanNumber_str + '.' + ParentCharge_str)
        spectrum_query.set('start_scan', ScanNumber_str)
        spectrum_query.set('end_scan', ScanNumber_str)
        spectrum_query.set('precursor_neutral_mass', psm_obj.measuredmass)
        spectrum_query.set('assumed_charge', ParentCharge_str)
        spectrum_query.set('index', str(iIndexUnique))
        iIndexUnique += 1
                    
        search_result = SubElement(spectrum_query, 'search_result')
        for idx, oPepScores in enumerate(psm_obj.list_list[op]):
            if idx == len(psm_obj.list_list[op]) - 1:
                break
            search_hit = SubElement(search_result, 'search_hit')
            search_hit.set('hit_rank', str(idx))
            search_hit.set('peptide', oPepScores.originalPep[1:-1])
            search_hit.set('protein', oPepScores.proteinlist[0])
            search_hit.set('num_tot_proteins', str(len(oPepScores.proteinlist)))
            search_hit.set('calc_neutral_pep_mass', oPepScores.calculatedmass)
            search_hit.set('massdiff', str(get_mass_diff(float(psm_obj.measuredmass), float(oPepScores.calculatedmass))))
            search_hit.set('num_tol_term', '2')
            search_hit.set('num_missed_cleavages', get_num_missed_cleavages(oPepScores.identifiedPep, cleave_residues_str))
            # alternative_protein
            if len(oPepScores.proteinlist) > 1:
                for iProteins in range(1, len(oPepScores.proteinlist)):
                    alternative_protein = SubElement(search_hit, 'alternative_protein')
                    alternative_protein.set('protein', oPepScores.proteinlist[iProteins])
            # modification_info
            local_modification_dict = get_modification_info(oPepScores.identifiedPep[1:-1], modification_dict)
            if local_modification_dict:
                modification_info = SubElement(search_hit, "modification_info")
                for key, value in local_modification_dict.iteritems():
                    mod_aminoacid_mass = SubElement(modification_info, 'mod_aminoacid_mass')
                    mod_aminoacid_mass.set('position', key)
                    mod_aminoacid_mass.set('mass', str(value))
            # search_score
            search_score = SubElement(search_hit, 'search_score')
            search_score.set('name', score_list_str[op])
            search_score.set('value', str(oPepScores.scorelist[op]))
            search_score = SubElement(search_hit, 'search_score')
            search_score.set('name', deltascore_str)
            search_score.set('value', str(oPepScores.scorediff[op]))

    # write into file
    document = ElementTree(root)
    document.write((filename +'_' +score_list_str[op]+ '.pep.xml'),
                   encoding='ISO-8859-1',
                   xml_declaration=True,
                   pretty_print=True)
    

def main(argv=None):
    
    # Display work start and time record
    start_time = datetime.now()
    sys.stderr.write('[%s] Beginning Sipros Ensemble Tabulating (%s)\n' % (curr_time(), get_version()))
    sys.stderr.write('------------------------------------------------------------------------------\n')

    if argv is None:
        argv = sys.argv
    # parse options
    (input_folder, _sConfig, output_folder, pepxml_output) = parse_options(argv)
    
    # multiple threading
    iNumThreads = cpu_count()
    # iNumThreads = 3
    if iNumThreads < 3:
        sys.stderr.write('This script needs at least three cores per CPU.')
    else:
        sys.stderr.write('%d threads will be used. \n' % iNumThreads)
    # Parse options and get config file
    sys.stderr.write('[Step 1] Parse options and get config file: Running -> ')
    # Call parse_config to open and read config file
    (config_dict, modification_dict, element_modification_list_dict) = parse_config(_sConfig)
    '''
    writePepxml('/media/xgo/Seagate/Proteomics/Experiments/Pepxml/OSU_D10_FASP_Elite_03202014_05.SE_Spe2Pep.txt.tab', config_dict, modification_dict, element_modification_list_dict)
    return
    '''
    # Get base_out for output
    base_out_default = 'Sipros_searches'
    sipros_file_list = get_file_list_with_ext(input_folder, 'Spe2Pep.txt')
    base_out = get_base_out(sipros_file_list, base_out_default, output_folder)
    
    iQueueSize = 10000
    qPsmUnprocessed = Queue(iQueueSize)
    qPsmProcessed = Queue(iQueueSize)
    sys.stderr.write('Done!\n')
    
    
    sys.stderr.write('[Step 2] Generate PSM table:                Running -> ')
    # File Reader (Producer)
    DataReader = Spe2PepReader(queue=qPsmUnprocessed, 
                               name='FileReader', 
                               searchname=config_dict[search_name_str],
                               inputFolder=input_folder)
    DataReader.daemon = True
    DataReader.start()
    
    # PSM processor (Consumer)
    for i in range(iNumThreads - 2):
        PsmProcessor = RankPsm(qPsmUnprocessed, qPsmProcessed, name=('PsmProcessor' + str(i)))
        PsmProcessor.daemon = True
        PsmProcessor.start()
    
    base_out_filename = base_out.split('/')[-1]
    base_out = output_folder + base_out_filename
    writePsm(base_out + '.tab', qPsmProcessed, iNumThreads - 2, pepxml_bool = pepxml_output)
    sys.stderr.write('Done!\n')
    
    
    # base_out_filename = base_out.split('/')[-1]
    # base_out = output_folder + base_out_filename
    sys.stderr.write('[Step 3] Merge Protein list:                Running -> ')
    refine_proteins(base_out+'.tab')
    sys.stderr.write('Done!\n')
    #writePin(base_out + '.tab', qPsmProcessed, iNumThreads - 2)
    #write_PepXML(output_folder, qPsmProcessed, iNumThreads - 2, config_dict)

    # Generate pep xml files
    if pepxml_output:
        sys.stderr.write('[Step 4] Generate Pepxml:                   Running -> ')
        
        sys.stderr.write('Done!\n')
    
    # Time record, calculate elapsed time, and display work end
    finish_time = datetime.now()
    duration = finish_time - start_time
    sys.stderr.write('------------------------------------------------------------------------------\n')
    sys.stderr.write('[%s] Ending Sipros Ensemble Tabulating\n' % curr_time())
    sys.stderr.write('Run complete [%s elapsed]\n' %  format_time(duration))
    
    print base_out + '.tab'

if __name__ == '__main__':
    sys.exit(main())
