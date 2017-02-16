#!/usr/bin/python

"""
sipros_post_module.py

sipros_post_module.py module includes common classes, definitions,
and funcitons for sipros post-processing programs.

Created by Tae-Hyuk (Ted) Ahn on 10/10/2012.
Copyright (c) 2012 Tae-Hyuk Ahn (ORNL). Allrights reserved.

Modified by Xuan Guo on 07/13/2016
"""

# # Import standard modules
import sys, os, re, math
from datetime import datetime
from collections import namedtuple
from multiprocessing import Process

# # Class Usage
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


# # Class for ignoring comments '#' in sipros file
class CommentedFile:
    def __init__(self, f, comment_string="#"):
        self.f = f
        self.comment_string = comment_string
    def next(self):
        line = self.f.next()
        while line.startswith(self.comment_string):
            line = self.f.next()
        return line
    def __iter__(self):
        return self


# # Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)


# # Returns the current time in a nice format
def curr_time():
    curr_time = datetime.now()
    return curr_time.strftime("%c")


# # Format time as a pretty string
def format_time(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)


# # Find string between two substrings
def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""


# # Division error handling
def divide(x, y):
    try:
        result = x / y
    except ZeroDivisionError as detail:
        print >> sys.stderr, 'Handling run-time error:', detail
        die('Program exit!')
    else:
        return result


# # Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as _f: pass
    except IOError as _e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")


# # Get file(s) list in working dir with specific file extension
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

        if len(file_list) == 0:
            print >> sys.stderr, "\nCannot open %s file(s)." % (file_ext)
            die("Program exit!")
        file_list = sorted(file_list)

    else:
        print >> sys.stderr, "\nCannot open working directory", working_dir
        die("Program exit!")

    return file_list

# # Get base output filename with input file list and base_out_default
def get_base_out(file_list, base_out_default, working_dir):

    # Get base output with common prefix
    base_out = os.path.commonprefix(file_list)
    base_out_filename = base_out.split('/')[-1]

    # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
    base_out = base_out.replace(".pep.txt", "_")

    # If base_out file name is less than 5, then use default baseout
    if len(base_out_filename) < 5:
        base_out = working_dir + base_out_default

    # If base common prefix ends with '_' or '.', then remove
    base_out = base_out[:-1] if (base_out[-1] in ('_', '.')) else base_out

    return base_out


# # list_to_string
# # if single element, then just convert to string
# # if multiple elements, then bracket {A,B}
def list_to_string(input_list):

    if len(input_list) > 1:
        converted_str = '{' + ','.join(input_list) + '}'
    else:
        converted_str = ''.join(input_list)

    return converted_str


# # list_to_bracket
# # bracket the list
def list_to_bracket(input_list):

    converted_str = '{' + ','.join(input_list) + '}'

    return converted_str


# # Class for sipros fields object
class SiprosFields(namedtuple('SiprosFields',
        ['Filename',
        'ScanNumber',
        'ParentCharge',
        'MeasuredParentMass',
        'CalculatedParentMass',
        'ScanType',
        'SearchName',
        'ScoringFunction',
        'Rank',
        'Score',
        'IdentifiedPeptide',
        'OriginalPeptide',
        'ProteinNames'])):
    def __init__(self):
        self.data = self

# # Class for spectrum fields object
class SpectrumFields(namedtuple('SpectrumFields',
        ['x',
        'Filename',
        'ScanNumber',
        'ParentCharge',
        'MeasuredParentMass',
        'ScanType',
        'SearchName',
        'TotalIntensity',
        'MaxIntensity'])):
    def __init__(self):
        self.data = self

# # Class for psm fields object
class PsmFields(namedtuple('PsmFields',
        ['x',
        'OriginalPeptide',
        'IdentifiedPeptide',
        'CalculatedParentMass',
        'MVH',
        'Xcorr',
        'WDP',
        'ProteinNames'])):
    def __init__(self):
        self.data = self


# # Class for sipros4 fields object
class Sipros4Fields(namedtuple('SiprosFields',
        ['Filename',
        'ScanNumber',
        'ParentCharge',
        'MeasuredParentMass',
        'CalculatedParentMass',
        'ScanType',
        'SearchName',
        'Rank',
        'MVH',
        'Xcorr',
        'WeightDotSum',
        'IdentifiedPeptide',
        'OriginalPeptide',
        'ProteinNames',
        'AveAtom',
        'StdAtom'])):

    def __init__(self):
        self.data = self


# # Class for PsmOutFields object
class PsmOutFields(namedtuple('PsmOutFields',
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',  # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',  # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',  # The difference score between the rank 1 and 2
         'DeltaP',  # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch'])):
    def __init__(self):
        self.data = self


# # Class for PsmOutFields object (for sipro4)
class Psm4OutFields(namedtuple('PsmOutFields',
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',  # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',  # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',  # The difference score between the rank 1 and 2
         'DeltaP',  # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch',
         'AveAtom',
         'StdAtom'])):
    def __init__(self):
        self.data = self


# # Class for PepSubFields object
class PepSubFields(namedtuple('PepSubFields',
        ['IdentifiedPeptide',
         'ParentCharge',
         'OriginalPeptide',
         'ProteinNames',
         'Score',
         'Filename',
         'ScanNumber',
         'ScanType',
         'SearchName'])):

    def __init__(self):
        self.data = self


# # Class for PepSubFields object
class PepDataFields(namedtuple('PepDataFields',
        ['IdentifiedPeptide',
         'ParentCharge',
         'BestScore',
         'ProteinNames'])):

    def __init__(self):
        self.data = self


# # Class for PepOutFields object
class PepOutFields(namedtuple('PepOutFields',
        ['IdentifiedPeptide',  # 0
         'ParentCharge',  # 1
         'OriginalPeptide',  # 2
         'ProteinNames',  # 3
         'ProteinCount',  # 4
         'TargetMatch',  # 5
         'SpectralCount',  # 6 number of PSMs matched to this peptide
         'BestScore',  # 7 the highest score of those PSMs
         'PSMs',  # 8 a list of PSMs matched to this peptide. Use {Filename[ScanNumber],Filename[ScanNumber]} format
         'ScanType',  # 9 ScanType
         'SearchName'])):  # 10 SearchName

    def __init__(self):
        self.data = self


# # Class for pretty float
class PrettyFloat(float):
    def __repr__(self):
        return "%0.5f" % self


# # A range function, that does accept float increments
def frange(start, end=None, inc=None):

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        fnext = start + len(L) * inc
        if inc > 0 and fnext >= end:
            break
        elif inc < 0 and fnext <= end:
            break
        L.append(fnext)

    return L


# # check sub list
def check_sub_list(list_A, list_B):

    check_status = True

    for list_A_item in list_A:
        if list_A_item not in list_B:
            check_status = False
        else:
            continue

    return check_status


# # get item list from parenthesis string as {AA,BB}
def get_item_list(input_string):

    input_string = input_string[1:-1]
    item_list = re.split(r"\s*[,]\s*", input_string.strip())

    return item_list


# # get_protein_count
def get_protein_count(protein_names):

    input_string = protein_names[1:-1]
    item_list = re.split(r"\s*[,]\s*", input_string.strip())
    protein_count = len(item_list)

    return protein_count


# # set float digit
def set_float_digit(input_val):

    if input_val is float:
        output_val = str("{0:.5f}".format(round(input_val, 5)))
    else:
        output_val = str(input_val)

    return output_val

# # peptide delete residues
def peptide_delete_residues(peptide_string):

    try:
        left_braket_index = peptide_string.index('[')
        right_braket_index = peptide_string.index(']')
        if len(peptide_string) > right_braket_index + 1:
            if peptide_string[right_braket_index + 1].isalpha():
                peptide_output = peptide_string[left_braket_index:right_braket_index + 1]
            else:
                peptide_output = peptide_string[left_braket_index:right_braket_index + 2]
        else:
            peptide_output = peptide_string[left_braket_index:right_braket_index + 1]

        return peptide_output
    except AttributeError:
        print >> sys.stderr, '\nCannot parse peptide correctly.\n'
        die("Program exit!")


# # merge protein names
def merge_protein_names(first_protein_names, second_protein_names):

    first_protein_list = get_item_list(first_protein_names)
    second_protein_list = get_item_list(second_protein_names)

    merge_protein_list = list(set(first_protein_list + second_protein_list))

    merge_protein_names = list_to_bracket(merge_protein_list)

    return merge_protein_names


# # Task wrapper
class PsmPack:

    def __init__(self, _iSize=1000, _iStartScanNumber=0):
        self.iSize = _iSize
        self.lPsm = []
        for _i in range(_iSize):
            self.lPsm.append([])
        self.iStartScanNumber = _iStartScanNumber
        self.bEmpty = True
        self.current = 0

    def add(self, lOnePsm, iScanNumber):
        self.lPsm[iScanNumber - self.iStartScanNumber].append(lOnePsm)
        self.bEmpty = False

    def empty(self):
        return self.bEmpty

    def __iter__(self):
        self.current = 0
        return self

    def next(self):
        if self.current >= self.iSize:
            raise StopIteration
        else:
            while not self.lPsm[self.current]:
                self.current += 1
                if self.current >= self.iSize:
                    raise StopIteration
            self.current += 1
            return self.lPsm[self.current - 1]


# # the mass difference, inverted, larger better.
def MassDiff(oPepScores):
    fDiff = abs(oPepScores.fCalculatedParentMass - oPepScores.fMeasuredParentMass)
    return -fDiff

# # the PTM score, the count of PTM, larger better
def PtmScore(oPepScores):
    s1 = ''.join([char if char.isalnum() else '$' for char in oPepScores.sIdentifiedPeptide ])
    return -(s1.count('$') - 2)

# # Rank Product
def RankProductInvert(liRank):
    fProduct = 1.0
    for i in liRank:
        fProduct *= i
    return 1 / fProduct

# # Pep info in the Spe2Pep file
class PepScores:

    def __init__(self, _fMeasuredParentMass, _iCharge, _sSearchName, sPeptideLine):
        self.fMeasuredParentMass = _fMeasuredParentMass
        self.iCharge = _iCharge
        asWords = sPeptideLine.split('\t')
        self.sIdentifiedPeptide = peptide_delete_residues(asWords[1])
        self.sOriginalPeptide = peptide_delete_residues(asWords[2])
        # self.sIdentifiedPeptide = peptide_delete_residues(asWords[2])
        # self.sOriginalPeptide = peptide_delete_residues(asWords[1])
        self.fCalculatedParentMass = float(asWords[3])
        self.lfScores = []
        self.liRanks = []
        self.lfScoreDiff = []  # difference between the current one with the second best one
        for e in asWords[4:-1]:
            self.lfScores.append(float(e))
        # remove the {}
        self.sProteinNames = (asWords[-1])[1:-1]
        self.fRankProduct = 0.0
        self.iRank = 0
        self.sSearchName = _sSearchName
        # delta -> comet way
        # diff -> difference between current one to the next best one
        # diffnor -> difference between current one to the next best one normalized by the current one
        self.lfDeltaRankProduct = []
        self.lfDeltaRankScore = []
        self.lfDiffRankProduct = []
        self.lfDiffRankScore = []
        self.lfDiffNorRankProduct = []
        self.lfDiffNorRankScore = []
        self.DeltaP = 'NA'
        if len(self.sOriginalPeptide) != len(self.sIdentifiedPeptide):
            self.DeltaP = 1

def numberTopRanks(liRanks):
    iCount = 0
    for i in liRanks:
        if i == 1:
            iCount += 1
    return iCount

def zero_divide(a, b):
    if b != 0:
        return a / b
    else:
        return 0

# # PSM info in the Spe2Pep file
class PepSpectrumMatch:

    iPurgeTrigger = 100
    iRankSave = 5

    def __init__(self, sSpectrumLine):
        asWords = sSpectrumLine.split('\t')
        self.sFileName = asWords[1]
        self.iScanNumber = int(asWords[2])
        self.iParentCharge = int(asWords[3])
        self.fMeasuredParentMass = float(asWords[4])
        self.sScanType = asWords[5]
        self.sSearchName = asWords[6]
        self.fTotalIntensity = float(asWords[7])
        self.fMaxIntensity = float(asWords[8])
        self.sRTime = '-1.000'
        if len(asWords) >= 10:
            self.sRTime = asWords[9]
        self.lPepScores = []
        self.oBestPep = None
        self.oSecondBestPep = None
        self.oRestPep = None
        self.lTopPep = []

    def addPepScores(self, pep):
        for e in self.lPepScores:
            if e.sIdentifiedPeptide == pep.sIdentifiedPeptide:
                words = pep.sProteinNames.split(',')
                for sProtein in words:
                    if e.sProteinNames.find(sProtein) == -1:
                        e.sProteinNames += ','
                        e.sProteinNames += sProtein
                return
        self.lPepScores.append(pep)
        if len(self.lPepScores) > self.iPurgeTrigger:
            self.purge()


    def purge(self):
        iNumScores = len(self.lPepScores[0].lfScores)
        for j in self.lPepScores:
            del j.liRanks[:]
        for i in range(iNumScores):
            lPep = sorted(self.lPepScores, key=lambda pep: (pep.lfScores[i], MassDiff(pep), PtmScore(pep), pep.sIdentifiedPeptide), reverse=True)
            iRank = 1
            for j in lPep:
                j.liRanks.append(iRank)
                iRank += 1
        liRanksNew = []
        for j in self.lPepScores:
            if any(i <= self.iRankSave for i in j.liRanks):
                liRanksNew.append(j)
        self.lPepScores = liRanksNew

    def ranking(self):
        del self.lTopPep[:]
        for j in self.lPepScores:
            del j.liRanks[:]
            del j.lfScoreDiff[:]
        iNumScores = len(self.lPepScores[0].lfScores)
        pep_rank_list = []
        for i in range(iNumScores):
            lPep = sorted(self.lPepScores, key=lambda pep: (pep.lfScores[i], MassDiff(pep), PtmScore(pep), pep.sIdentifiedPeptide), reverse=True)
            pep_rank_list.append(lPep)
            iRank = 1
            for j in lPep:
                j.liRanks.append(iRank)
                iRank += 1
            # score rank -> score differential
            for j in range(0, len(lPep) - 1):
                lPep[j].lfDeltaRankScore.append(1 - zero_divide(lPep[j + 1].lfScores[i], lPep[0].lfScores[i]))
                '''
                if j == 0:
                    lPep[j].lfDeltaRankScore.append(1.0 - zero_divide(lPep[1].lfScores[i], lPep[0].lfScores[i]))
                else:
                    lPep[j].lfDeltaRankScore.append(zero_divide(lPep[j].lfScores[i], lPep[0].lfScores[i]) - 1.0)
                '''
                lPep[j].lfDiffRankScore.append(lPep[j].lfScores[i] - lPep[j + 1].lfScores[i])
                lPep[j].lfDiffNorRankScore.append(zero_divide(lPep[j].lfScores[i] - lPep[j + 1].lfScores[i], lPep[j].lfScores[i]))
            lPep[len(lPep) - 1].lfDeltaRankScore.append(0)
            # lPep[len(lPep) - 1].lfDeltaRankScore.append(zero_divide(lPep[len(lPep) - 1].lfScores[i], lPep[0].lfScores[i]) - 1.0)
            lPep[len(lPep) - 1].lfDiffRankScore.append(0)
            lPep[len(lPep) - 1].lfDiffNorRankScore.append(0)
            
            self.lTopPep.append(lPep[0])
            if len(lPep) >= 2:
                for j in lPep:
                    j.lfScoreDiff.append(j.lfScores[i] - lPep[1].lfScores[i])
            else:
                j.lfScoreDiff.append(0.0)

        # rank product -> score differential
        for i in self.lPepScores:
            i.fRankProduct = RankProductInvert(i.liRanks)
        # lPep: ranked by rank product
        lPep = sorted(self.lPepScores, key=lambda pep: (pep.fRankProduct, MassDiff(pep), PtmScore(pep), pep.sIdentifiedPeptide), reverse=True)
        for j in range(0, len(lPep) - 1):
            lPep[j].iRank = j
            for i in range(iNumScores):
                lPep[j].lfDeltaRankProduct.append(1 - zero_divide(lPep[j + 1].lfScores[i], lPep[0].lfScores[i]))
                lPep[j].lfDiffRankProduct.append(lPep[j].lfScores[i] - lPep[j + 1].lfScores[i])
                lPep[j].lfDiffNorRankProduct.append(zero_divide(lPep[j].lfScores[i] - lPep[j + 1].lfScores[i], lPep[j].lfScores[i]))
        lPep[-1].iRank = len(lPep) - 1
        
        for i in range(iNumScores):
            lPep[len(lPep) - 1].lfDeltaRankProduct.append(0)
            lPep[len(lPep) - 1].lfDiffRankProduct.append(0)
            lPep[len(lPep) - 1].lfDiffNorRankProduct.append(0)
        
        for i in self.lTopPep:
            if numberTopRanks(i.liRanks) >= 2:
                self.oBestPep = i
                self.iParentCharge = self.oBestPep.iCharge
                self.sSearchName = self.oBestPep.sSearchName
                if len(lPep) > 1:
                    if self.oBestPep == lPep[0]:
                        self.oSecondBestPep = lPep[1]
                    else:
                        self.oSecondBestPep = lPep[0]
                    # PTM in the best pep
                    if len(i.sIdentifiedPeptide) != len(i.sOriginalPeptide):
                        for pep in lPep:
                            if pep.sIdentifiedPeptide != i.sIdentifiedPeptide and \
                            abs(pep.fCalculatedParentMass - i.fCalculatedParentMass) < 0.00005 and \
                            abs(math.fabs(pep.fMeasuredParentMass - pep.fCalculatedParentMass) - math.fabs(i.fMeasuredParentMass - i.fCalculatedParentMass)) < 0.00005 and \
                            pep.sOriginalPeptide == i.sOriginalPeptide and \
                            len(pep.sIdentifiedPeptide) == len(i.sIdentifiedPeptide):
                                avg_deltaP = 0.0
                                for s in range(iNumScores):
                                    if self.lTopPep[s].lfScores[s] != 0:
                                        avg_deltaP += (i.lfScores[s] - pep.lfScores[s])/self.lTopPep[s].lfScores[s]
                                avg_deltaP /= float(iNumScores)
                                self.oBestPep.DeltaP = avg_deltaP
                                return
                return
        # if SA == 1, calculate DeltaP individually
        for s in range(iNumScores):
            for lPep_local in pep_rank_list:
                # contain PTM
                if len(lPep_local[0].sIdentifiedPeptide) != len(lPep_local[0].sOriginalPeptide):
                    for pep in lPep_local:
                        if pep.sIdentifiedPeptide != lPep_local[0].sIdentifiedPeptide and \
                        abs(pep.fCalculatedParentMass - lPep_local[0].fCalculatedParentMass) < 0.00005 and \
                        abs(math.fabs(pep.fMeasuredParentMass - pep.fCalculatedParentMass) - math.fabs(lPep_local[0].fMeasuredParentMass - lPep[0].fCalculatedParentMass)) < 0.00005 and \
                        pep.sOriginalPeptide == lPep_local[0].sOriginalPeptide and \
                        len(pep.sIdentifiedPeptide) == len(lPep_local[0].sIdentifiedPeptide):
                            avg_deltaP = 0.0
                            for si in range(iNumScores):
                                if self.lTopPep[si].lfScores[si] != 0:
                                    avg_deltaP += (lPep_local[0].lfScores[si] - pep.lfScores[si])/self.lTopPep[si].lfScores[si]
                            avg_deltaP /= float(iNumScores)
                            lPep_local[0].DeltaP = avg_deltaP
                            break
        self.oBestPep = lPep[0]
        self.iParentCharge = self.oBestPep.iCharge
        self.sSearchName = self.oBestPep.sSearchName
        if len(lPep) > 1:
            if self.oBestPep == lPep[0]:
                self.oSecondBestPep = lPep[1]
            else:
                self.oSecondBestPep = lPep[0]
        # debug
        # if protein_type(self.oBestPep.sProteinNames) == -1:
            # pass

    def __repr__(self):
        # Filename    ScanNumber    ParentCharge    MeasuredParentMass    ScanType    SearchName
        # IdentifiedPeptide    OriginalPeptide    CalculatedParentMass     MVH    WDP    Xcorr    ProteinNames
        return "%s\t%d\t%d\t%f\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t{%s}" % (self.sFileName,
                                                                       self.iScanNumber,
                                                                       self.iParentCharge,
                                                                       self.fMeasuredParentMass,
                                                                       self.sScanType,
                                                                       self.sSearchName,
                                                                       self.oBestPep.sIdentifiedPeptide,
                                                                       self.oBestPep.sOriginalPeptide,
                                                                       self.oBestPep.fCalculatedParentMass,
                                                                       self.oBestPep.lfScores[0],
                                                                       self.oBestPep.lfScores[1],
                                                                       self.oBestPep.lfScores[2],
                                                                       self.oBestPep.sProteinNames)

    def get_other_features(self):
        # score differential
        # comet way
        # difference
        # difference normalization
        if len(self.oBestPep.lfDeltaRankProduct) == 0:
            return '0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000'
        
        return "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (self.oBestPep.lfDeltaRankProduct[0],
                                                                                self.oBestPep.lfDeltaRankProduct[1],
                                                                                self.oBestPep.lfDeltaRankProduct[2],
                                                                                self.oBestPep.lfDeltaRankScore[0],
                                                                                self.oBestPep.lfDeltaRankScore[1],
                                                                                self.oBestPep.lfDeltaRankScore[2],
                                                                                self.oBestPep.lfDiffRankProduct[0],
                                                                                self.oBestPep.lfDiffRankProduct[1],
                                                                                self.oBestPep.lfDiffRankProduct[2],
                                                                                self.oBestPep.lfDiffRankScore[0],
                                                                                self.oBestPep.lfDiffRankScore[1],
                                                                                self.oBestPep.lfDiffRankScore[2],
                                                                                self.oBestPep.lfDiffNorRankProduct[0],
                                                                                self.oBestPep.lfDiffNorRankProduct[1],
                                                                                self.oBestPep.lfDiffNorRankProduct[2],
                                                                                self.oBestPep.lfDiffNorRankScore[0],
                                                                                self.oBestPep.lfDiffNorRankScore[1],
                                                                                self.oBestPep.lfDiffNorRankScore[2])
        
    def get_second_pep(self):
        return "%s\t%d\t%d\t%f\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t{%s}" % (self.sFileName,
                                                                       self.iScanNumber,
                                                                       self.oSecondBestPep.iCharge,
                                                                       self.fMeasuredParentMass,
                                                                       self.sScanType,
                                                                       self.oSecondBestPep.sSearchName,
                                                                       self.oSecondBestPep.sIdentifiedPeptide,
                                                                       self.oSecondBestPep.sOriginalPeptide,
                                                                       self.oSecondBestPep.fCalculatedParentMass,
                                                                       self.oSecondBestPep.lfScores[0],
                                                                       self.oSecondBestPep.lfScores[1],
                                                                       self.oSecondBestPep.lfScores[2],
                                                                       self.oSecondBestPep.sProteinNames)
    
    def get_other_features_second_pep(self):
        if len(self.oSecondBestPep.lfDeltaRankProduct) == 0:
            return '0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000\t0.000'
        return "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (self.oBestPep.lfDeltaRankProduct[0],
                                                                                self.oSecondBestPep.lfDeltaRankProduct[1],
                                                                                self.oSecondBestPep.lfDeltaRankProduct[2],
                                                                                self.oSecondBestPep.lfDeltaRankScore[0],
                                                                                self.oSecondBestPep.lfDeltaRankScore[1],
                                                                                self.oSecondBestPep.lfDeltaRankScore[2],
                                                                                self.oSecondBestPep.lfDiffRankProduct[0],
                                                                                self.oSecondBestPep.lfDiffRankProduct[1],
                                                                                self.oSecondBestPep.lfDiffRankProduct[2],
                                                                                self.oSecondBestPep.lfDiffRankScore[0],
                                                                                self.oSecondBestPep.lfDiffRankScore[1],
                                                                                self.oSecondBestPep.lfDiffRankScore[2],
                                                                                self.oSecondBestPep.lfDiffNorRankProduct[0],
                                                                                self.oSecondBestPep.lfDiffNorRankProduct[1],
                                                                                self.oSecondBestPep.lfDiffNorRankProduct[2],
                                                                                self.oSecondBestPep.lfDiffNorRankScore[0],
                                                                                self.oSecondBestPep.lfDiffNorRankScore[1],
                                                                                self.oSecondBestPep.lfDiffNorRankScore[2])
    
    def all_top_ranked_psm(self):
        str_list = []
        for pep in self.lTopPep:
            feature_list = []
            feature_list.append(self.sFileName)
            feature_list.append(str(self.iScanNumber))
            feature_list.append(str(self.iParentCharge))
            feature_list.append(str(self.fMeasuredParentMass))
            feature_list.append(self.sScanType)
            feature_list.append(self.sSearchName)
            feature_list.append(pep.sIdentifiedPeptide)
            feature_list.append(pep.sOriginalPeptide)
            feature_list.append(str(pep.fCalculatedParentMass))
            feature_list.extend((str(x) for x in pep.lfScores))
            feature_list.append('{'+pep.sProteinNames+'}')
            feature_list.append(str(num_agreement(pep.liRanks)))
            feature_list.extend((str(x) for x in pep.lfDeltaRankProduct))
            feature_list.extend((str(x) for x in pep.lfDeltaRankScore))
            feature_list.extend((str(x) for x in pep.lfDiffRankProduct))
            feature_list.extend((str(x) for x in pep.lfDiffRankScore))
            feature_list.extend((str(x) for x in pep.lfDiffNorRankProduct))
            feature_list.extend((str(x) for x in pep.lfDiffNorRankScore))
            feature_list.append(self.sRTime)
            feature_list.append((str(pep.iRank)))
            feature_list.append((str(pep.DeltaP)))
            str_list.append('\t'.join(feature_list))
        return '\n'.join(str_list)
    
    def pin(self):
        lPin = []
        '''
        id    label    ScanNr  
        ScoreAgreement  
        deltCn    Xcorr    Mass    
        PepLen    Charge1    Charge2    
        Charge3    Charge4    Charge5    
        Charge6    enzInt    dM    absdM    
        peptide    proteinId1
        '''
        lPin.append((self.sFileName
                     + '_'
                     + str(self.iScanNumber)
                     + '_'
                     + str(self.iParentCharge)))
        lProteins = []
        iType = protein_type(self.oBestPep.sProteinNames, lProteins)
        if iType == 2:
            lPin.append('-1')
        else:
            lPin.append('1')
            lProteins = self.removeReverse(lProteins)
        lPin.append(str(self.iScanNumber))
        lPin.append(str(num_agreement(self.oBestPep.liRanks)))
        lPin.append(str(self.oBestPep.lfScoreDiff[0]))
        lPin.append(str(self.oBestPep.lfScores[0]))
        lPin.append(str(math.log(self.oBestPep.liRanks[0])))
        lPin.append(str(self.oBestPep.lfScoreDiff[1]))
        lPin.append(str(self.oBestPep.lfScores[1]))
        lPin.append(str(math.log(self.oBestPep.liRanks[1])))
        lPin.append(str(self.oBestPep.lfScoreDiff[2]))
        lPin.append(str(self.oBestPep.lfScores[2]))
        lPin.append(str(math.log(self.oBestPep.liRanks[2])))
        lPin.append(str(self.fMeasuredParentMass))
        lPin.append(str(len(self.oBestPep.sOriginalPeptide) - 2))
        for i in range(1, 6, 1):
            if self.iParentCharge == i:
                lPin.append('1')
            else:
                lPin.append('0')
        if self.iParentCharge >= 6:
            lPin.append('1')
        else:
            lPin.append('0')
        lPin.append(str(self.enzInt()))
        fTemp = self.fMeasuredParentMass - self.oBestPep.fCalculatedParentMass
        lPin.append(str(fTemp))
        lPin.append(str(abs(fTemp)))
        lPin.append(self.percolatorPeptide())
        # lPin.extend(lProteins)
        lPin.append(lProteins[0])
        return ('\t'.join(lPin), len(lProteins))

    def removeReverse(self, lProteins):
        newlist = []
        for sProtein in lProteins:
            if not sProtein.startswith('Rev_'):
                newlist.append(sProtein)
        return newlist

    def percolatorPeptide(self):
        return 'K.' + self.oBestPep.sIdentifiedPeptide[1:-1] + '.A'

    def enzInt(self):
        iNumIntEnz = -1
        iNumIntEnz += self.oBestPep.sIdentifiedPeptide.count('K')
        iNumIntEnz += self.oBestPep.sIdentifiedPeptide.count('R')
        return iNumIntEnz

# # lOnePsm: + spectrum * peptide * peptide + spectrum
def SelectTopRankedPsm(lOnePsm):
    psm = PepSpectrumMatch(lOnePsm[0][0])
    '''
    if psm.iScanNumber == 18472:
        print 'check'
    '''
    iCharge = 0
    sSearchName = ''
    for PsmInOneFile in lOnePsm:
        for sline in PsmInOneFile:
            if sline[0] == '+':
                # a spectrum line
                iCharge = get_charge(sline)
                sSearchName = get_search_name(sline)
            else:
                # a peptide line
                pep = PepScores(psm.fMeasuredParentMass, iCharge, sSearchName, sline)
                psm.addPepScores(pep)
    # sorting and then ranking
    psm.ranking()
    return psm

# # peak a line from the file
def peek_line(f):
    pos = f.tell()
    sline = f.readline()
    f.seek(pos)
    return sline

# # get the scan number from a line
def get_scan_number(sLine, sDelimiter='\t', iFirstDelimiter=2):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    iScanNumber = int(sLine[iBegin:iPos])
    return iScanNumber

# # get the charge from a line
def get_charge(sLine, sDelimiter='\t', iFirstDelimiter=3):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    iCharge = int(sLine[iBegin:iPos])
    return iCharge

# # get the search name from a line
def get_search_name(sLine, sDelimiter='\t', iFirstDelimiter=6):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    sSearchName = sLine[iBegin:iPos]
    return sSearchName

# # get the PSM with scan number less than the upper scan number bound
def get_psm(f, lPsm, sSpectrum='+', sPeptide='*', iUpperScanNumber=0):
    bEof = True
    lOnePsm = []
    while True:
        pos = f.tell()
        sline = f.readline().strip()
        # # end of file
        if not sline:
            break
        iScanNumber = 0
        if sline[0] == sSpectrum:
            iScanNumber = get_scan_number(sline)
            if iScanNumber < iUpperScanNumber:
                bEof = False
                lOnePsm = []
                lOnePsm.append(sline)
                lPsm.add(lOnePsm, iScanNumber)
            else:
                # roll back the previous position
                f.seek(pos)
                bEof = False
                break
        else:
            lOnePsm.append(sline)
    return bEof

# # skip the comment area and the header
def skip_comment(f, sComment='#', iLineHeader=0):
    pos = f.tell()
    sline = f.readline()
    while(sline[0] == sComment):
        pos = f.tell()
        sline = f.readline()
    f.seek(pos)
    for _i in range(iLineHeader):
        f.readline()
    # print f.readline()

# # Spe2Pep reader
class Spe2PepReader(Process):

    def __init__(self, queue=None, name=None, searchname=None, inputFolder=None):
        super(Spe2PepReader, self).__init__()
        self.name = name
        self.qPsmUnprocessed = queue;
        self.iNumScanProcessed = 0
        self.sSearchName = searchname
        self.FileList = []
        self.iScanInterval = 1000
        self.sInputFolder = inputFolder

    # # list files with 'Spe2Pep.txt' extensions
    # # put the search results for the same FT2 file into a list
    def categorizeSpe2PepFile(self, sWorkingDirectory):
        lFileList = get_file_list_with_ext(sWorkingDirectory, 'Spe2Pep.txt')
        sFt2Name = ''
        lFt2Name = []
        iIndexFt2 = 0
        for sFileName in lFileList:
            iPos = sFileName.rfind(self.sSearchName)
            if iPos != -1:
                # a '.' is before the search name, so iPos-1
                sFt2Name = sFileName[0:iPos - 1]
                if sFt2Name in lFt2Name:
                    iIndexFt2 = lFt2Name.index(sFt2Name)
                    self.FileList[iIndexFt2].append(sFileName)
                else:
                    lFt2Name.append(sFt2Name)
                    lNewFt2 = []
                    lNewFt2.append(sFileName)
                    self.FileList.append(lNewFt2)

    def readData(self):
        for _id, lFiles in enumerate(self.FileList):
            lFileReader = []
            for sFiles in lFiles:
                oFile = open(sFiles, 'r')
                skip_comment(oFile, iLineHeader=2)
                lFileReader.append(oFile)
            # # peek the first scan number
            iSmallestScanNumer = sys.maxint
            for f in lFileReader:
                sLine = peek_line(f)
                iScanNumber = get_scan_number(sLine)
                if iScanNumber < iSmallestScanNumer:
                    iSmallestScanNumer = iScanNumber
            # #
            iLastScanExcluded = iSmallestScanNumer
            bReachEof = False
            while(not bReachEof):
                psmPack = PsmPack(_iSize=self.iScanInterval, _iStartScanNumber=iLastScanExcluded)
                iLastScanExcluded = iLastScanExcluded + self.iScanInterval
                bReachEof = True
                for f in lFileReader:
                    bReachEof = bReachEof & (get_psm(f, psmPack, iUpperScanNumber=iLastScanExcluded))
                if not psmPack.empty():
                    # # add to the job queue
                    for psm in psmPack:
                        self.qPsmUnprocessed.put(psm, True)
                        self.iNumScanProcessed += 1
                        if self.iNumScanProcessed % 100 == 0:
                            pass
                            # print 'Read # scans %d' % self.iNumScanProcessed
                            # print 'Queue size %f' % self.qPsmUnprocessed.qsize()

            # # close the file reader
            for f in lFileReader:
                f.close()

    def run(self):
        if not self.qPsmUnprocessed:
            sys.stderr.write("Job Queue error in PSM reading.")
        self.categorizeSpe2PepFile(self.sInputFolder)
        self.readData()
        self.qPsmUnprocessed.put(None)

# # thread class for ranking the PSM
class RankPsm(Process):

    def __init__(self, qPsmUnprocessed, qPsmProcessed, name=None):
        super(RankPsm, self).__init__()
        self.name = name
        self.qPsmUnprocessed = qPsmUnprocessed
        self.qPsmProcessed = qPsmProcessed
        self.iCount = 0
        return

    def run(self):
        if not self.qPsmUnprocessed:
            sys.stderr.write("Job Queue error in PSM ranking.")
        if not self.qPsmProcessed:
            sys.stderr.write("Job Queue error in PSM ranking.")
        while True:
            psm = self.qPsmUnprocessed.get(True)
            if psm is None:
                break
            oPsm = SelectTopRankedPsm(psm)
            del psm
            self.iCount += 1
            if self.iCount % 10 == 0:
                pass
                # print "Rank # scans %i" % self.iCount
            self.qPsmProcessed.put(oPsm, True)
        self.qPsmProcessed.put(None)
        self.qPsmUnprocessed.put(None)
        return


# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    asProteins = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        lProtein.extend(asProteins[:])
    for sProtein in asProteins:
        if not (sProtein.startswith('Rev_') or sProtein.startswith('Dec_')):
            return 1
    for sProtein in asProteins:
        if sProtein.startswith('Dec_'):
            return 3
    return 2

# # score agreement recode
def agreement(liRank):
    iIndex = 0
    for i in liRank:
        if i == 1:
            iIndex = (iIndex << 1) + 1;
        else:
            iIndex = (iIndex << 1)
    return iIndex

# # number of score agreement
def num_agreement(liRanks):
    iNum = 0
    for i in liRanks:
        if i == 1:
            iNum += 1
    return iNum

# # write the PSM table
def writePin(sOutputFile, qPsmProcessed, iNumRankers):
    iNumTarget = 0
    iNumReverse = 0
    iNumShuffle = 0
    iNumProcessedScans = 0
    liAgreeRecord = []
    liTarget = []
    liReverse = []
    liShuffle = []
    for _i in range(8):
        liAgreeRecord.append(0)
        liTarget.append(0)
        liReverse.append(0)
        liShuffle.append(0)
    iMaxNumProtein = 0
    with open(sOutputFile, 'w') as f:
        # PSM
        while True:
            psm = qPsmProcessed.get(True)
            if psm is None:
                iNumRankers -= 1
                if iNumRankers == 0:
                    break
                else:
                    continue
            (sLine, iNumProteins) = psm.pin()
            if iNumProteins > iMaxNumProtein:
                iMaxNumProtein = iNumProteins
            f.write(sLine)
            f.write('\n')
            iTemp2 = agreement(psm.oBestPep.liRanks)
            liAgreeRecord[iTemp2] += 1
            iTemp = protein_type(psm.oBestPep.sProteinNames)
            if iTemp == 1:
                iNumTarget += 1
                liTarget[iTemp2] += 1
            elif iTemp == 2:
                iNumReverse += 1
                liReverse[iTemp2] += 1
            else:
                iNumShuffle += 1
                liShuffle[iTemp2] += 1
            del psm
            iNumProcessedScans += 1
            if iNumProcessedScans % 100 == 0:
                print " Processed & Saved %i Scans\r" % iNumProcessedScans

    print "\nTarget #: %d\tReverse #: %d\tShuffle #: %d" % (iNumTarget, iNumReverse, iNumShuffle)
    print liAgreeRecord
    print liTarget
    print liReverse
    print liShuffle
    print 'Max number of proteins: %i' % iMaxNumProtein

bAdditionalFeatures = True
header_message_additional_feature = '\tDeltaRP1\tDeltaRP2\tDeltaRP3\tDeltaRS1\tDeltaRS2\tDeltaRS3\tDiffRP1\tDiffRP2\tDiffRP3\tDiffRS1\tDiffRS2\tDiffRS3\tDiffNorRP1\tDiffNorRP2\tDiffNorRP3\tDiffNorRS1\tDiffNorRS2\tDiffNorRS3'
bExtraNegativePsm = False
bRTime = True
bAdditionPepScoreAgreementNotThree = True

# # write the PSM table
def writePsm(sOutputFile, qPsmProcessed, iNumRankers):
    iNumTarget = 0
    iNumReverse = 0
    iNumShuffle = 0
    iNumProcessedScans = 0
    liAgreeRecord = []
    liTarget = []
    liReverse = []
    liShuffle = []
    for _i in range(8):
        liAgreeRecord.append(0)
        liTarget.append(0)
        liReverse.append(0)
        liShuffle.append(0)
    if bExtraNegativePsm:
        sOutputFileExtraNegativePsm = sOutputFile[:-4] + '_nagative.tab'
        fExtra = open(sOutputFileExtraNegativePsm, 'w')
    
    with open(sOutputFile, 'w') as f:
        # header
        f.write('FileName\t')
        f.write('ScanNumber\t')
        f.write('ParentCharge\t')
        f.write('MeasuredParentMass\t')
        f.write('ScanType\t')
        f.write('SearchName\t')
        f.write('IdentifiedPeptide\t')
        f.write('OriginalPeptide\t')
        f.write('CalculatedParentMass\t')
        f.write('MVH\t')
        f.write('Xcorr\t')
        f.write('WDP\t')
        f.write('ProteinNames\t')
        f.write('ScoreAgreement')
        if bAdditionalFeatures:
            f.write(header_message_additional_feature)
        if bRTime:
            f.write('\tRetentionTime')
        f.write('\n')
        
        if bExtraNegativePsm:
            fExtra.write('FileName\t')
            fExtra.write('ScanNumber\t')
            fExtra.write('ParentCharge\t')
            fExtra.write('MeasuredParentMass\t')
            fExtra.write('ScanType\t')
            fExtra.write('SearchName\t')
            fExtra.write('IdentifiedPeptide\t')
            fExtra.write('OriginalPeptide\t')
            fExtra.write('CalculatedParentMass\t')
            fExtra.write('MVH\t')
            fExtra.write('Xcorr\t')
            fExtra.write('WDP\t')
            fExtra.write('ProteinNames\t')
            fExtra.write('ScoreAgreement')
            if bAdditionalFeatures:
                fExtra.write(header_message_additional_feature)
            fExtra.write('\n')
        
        # PSM
        while True:
            psm = qPsmProcessed.get(True)
            if psm is None:
                iNumRankers -= 1
                if iNumRankers == 0:
                    break
                else:
                    continue
            
            if bAdditionPepScoreAgreementNotThree and num_agreement(psm.oBestPep.liRanks) < 2:
                f.write(psm.all_top_ranked_psm())
            else:
                f.write(repr(psm))
                f.write('\t')
            #f.write(str(agreement(psm.oBestPep.liRanks)))
                if bAdditionalFeatures:
                    f.write(str(num_agreement(psm.oBestPep.liRanks)))
                    f.write('\t')
                    f.write(psm.get_other_features())
                else:
                    f.write(str(num_agreement(psm.oBestPep.liRanks)))
                if bRTime:
                    f.write('\t')
                    f.write(psm.sRTime)
                f.write('\t')
                f.write(str(psm.oBestPep.iRank))
                f.write('\t')
                f.write(str(psm.oBestPep.DeltaP))
            
            f.write('\n')
            iTemp2 = agreement(psm.oBestPep.liRanks)
            liAgreeRecord[iTemp2] += 1
            iTemp = protein_type(psm.oBestPep.sProteinNames)
            if iTemp == 1:
                iNumTarget += 1
                liTarget[iTemp2] += 1
            elif iTemp == 2:
                iNumReverse += 1
                liReverse[iTemp2] += 1
            else:
                iNumShuffle += 1
                liShuffle[iTemp2] += 1
                
            if bExtraNegativePsm and psm.oSecondBestPep != None:
                fExtra.write(psm.get_second_pep())
                fExtra.write('\t')
                if bAdditionalFeatures:
                    fExtra.write(str(num_agreement(psm.oSecondBestPep.liRanks)))
                    fExtra.write('\t')
                    fExtra.write(psm.get_other_features())
                else:
                    fExtra.write(str(num_agreement(psm.oSecondBestPep.liRanks)))
                fExtra.write('\n')
                
            del psm
            iNumProcessedScans += 1
            if iNumProcessedScans % 100 == 0:
                print " Processed & Saved %i Scans\r" % iNumProcessedScans
    
    if bExtraNegativePsm:
        fExtra.close()
    
    print "\nTarget #: %d\tReverse #: %d\tShuffle #: %d" % (iNumTarget, iNumReverse, iNumShuffle)
    print liAgreeRecord
    print liTarget
    print liReverse
    print liShuffle

from lxml.etree import ElementTree, Element, SubElement

pep_iden_str = '[Peptide_Identification]'
search_name_str = 'Search_Name'
FASTA_Database_str = 'FASTA_Database'
Maximum_Missed_Cleavages_str = 'Maximum_Missed_Cleavages'
Cleave_After_Residues_str = 'Cleave_After_Residues'

def get_num_missed_cleavages(peptide_str, cleave_residues_str):
    num = 0
    for c in cleave_residues_str:
        num += peptide_str.count(c)

    return str(num - 1)

def get_modification_info(peptide_str, modification_label_dict):
    modification_dict = {}
    for key, value in modification_label_dict.iteritems():
        beg = -1
        beg = peptide_str.find(key, beg + 1)
        while beg != -1:
            modification_dict[str(beg - len(modification_dict))] = value
            beg = peptide_str.find(key, beg + 1)
    return modification_dict

def write_PepXML(output_folder, qPsmProcessed, iNumRankers, config_dict):

    input_filename = ''
    root = Element('msms_pipeline_analysis')
    iIndexUnique = 1
    iScoreId = 1  # 0: MVH, 1: Xcorr, 2: WDP
    score_list_str = ['mvh', 'xcorr', 'mvh']
    search_engine_str = ['MyriMatch', 'Comet', 'MyriMatch']
    iNumProcessedScans = 0
    cleave_residues_str = config_dict[Cleave_After_Residues_str]

    modification_label_dict = {'~':'15.994915', '!': '0.984016'}

    while True:
        psm = qPsmProcessed.get(True)
        if psm is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
        if input_filename != psm.sFileName.split('.')[0]:
            if input_filename != '':
                document = ElementTree(root)
                document.write((output_folder + input_filename + '.pep.xml'),
                               encoding='ISO-8859-1',
                               xml_declaration=True,
                               pretty_print=True)
                iIndexUnique = 1
            input_filename = psm.sFileName.split('.')[0]
            xmlns = "http://regis-web.systemsbiology.net/pepXML"
            xsi = "http://www.w3.org/2001/XMLSchema-instance"
            schemaLocation = "http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd"
            root = Element("{" + xmlns + "}msms_pipeline_analysis",
                           nsmap={'xsi':xsi, None:xmlns},
                           attrib={"{" + xsi + "}schemaLocation"  : schemaLocation})
            root.set('date', datetime.now().isoformat())
            root.set('summary_xml', (input_filename + '.pepXML'))
            # root.set('xmlns', "http://regis-web.systemsbiology.net/pepXML")
            # root.set('xmlns:xsi', "http://www.w3.org/2001/XMLSchema-instance")
            # root.set('xsi:schemaLocation', "http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd")

            if iScoreId != 1:  # Xcorr does not need analysis summary
                analysis_summary = SubElement(root, 'analysis_summary')
                analysis_summary.set('analysis', "Sipros")
                analysis_summary.set('version', "10")
                analysis_summary.set('time', (datetime.now().isoformat()))

            msms_run_summary = SubElement(root, 'msms_run_summary')
            msms_run_summary.set('base_name', input_filename)
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
            search_summary.set('base_name', input_filename)
            search_summary.set('search_engine', search_engine_str[iScoreId])
            search_summary.set('precursor_mass_type', 'monoisotopic')
            search_summary.set('fragment_mass_type', 'monoisotopic')
            search_summary.set('out_data_type', '')
            search_summary.set('out_data', '')

            search_database = SubElement(search_summary, 'search_database')
            search_database.set('local_path', config_dict[FASTA_Database_str])
            search_database.set('database_name', 'SDB')
            search_database.set('type', 'AA')

            enzymatic_search_constraint = SubElement(search_summary, 'enzymatic_search_constraint')
            enzymatic_search_constraint.set('enzyme', 'Trypsin/P')
            enzymatic_search_constraint.set('max_num_internal_cleavages', config_dict[Maximum_Missed_Cleavages_str])
            enzymatic_search_constraint.set('min_number_termini', '2')

            aminoacid_modification = SubElement(search_summary, 'aminoacid_modification')
            aminoacid_modification.set('aminoacid', 'M')
            aminoacid_modification.set('massdiff', '15.994915')
            aminoacid_modification.set('mass', '147.0353846062')
            aminoacid_modification.set('variable', 'Y')
            # aminoacid_modification.set('symbol', '~')

            aminoacid_modification = SubElement(search_summary, 'aminoacid_modification')
            aminoacid_modification.set('aminoacid', 'C')
            aminoacid_modification.set('massdiff', '57.021464')
            aminoacid_modification.set('mass', '160.030649')
            aminoacid_modification.set('variable', 'N')

            # aminoacid_modification = SubElement(search_summary, 'aminoacid_modification')
            # aminoacid_modification.set('aminoacid', 'N')
            # aminoacid_modification.set('massdiff', '0.984016')
            # aminoacid_modification.set('variable', 'Y')
            # aminoacid_modification.set('symbol', '!')

            # aminoacid_modification = SubElement(search_summary, 'aminoacid_modification')
            # aminoacid_modification.set('aminoacid', 'Q')
            # aminoacid_modification.set('massdiff', '0.984016')
            # aminoacid_modification.set('variable', 'Y')
            # aminoacid_modification.set('symbol', '!')

        # query results
        spectrum_query = SubElement(msms_run_summary, 'spectrum_query')
        ScanNumber_str = str(psm.iScanNumber)
        ParentCharge_str = str(psm.iParentCharge)
        spectrum_query.set('spectrum', input_filename + '.' + ScanNumber_str + '.' + ScanNumber_str + '.' + ParentCharge_str)
        spectrum_query.set('start_scan', ScanNumber_str)
        spectrum_query.set('end_scan', ScanNumber_str)
        spectrum_query.set('precursor_neutral_mass', str(psm.fMeasuredParentMass))
        spectrum_query.set('assumed_charge', ParentCharge_str)
        spectrum_query.set('index', str(iIndexUnique))

        search_result = SubElement(spectrum_query, 'search_result')
        for oPepScores in psm.lPepScores:
            if oPepScores.liRanks[iScoreId] > 5:
                continue
            search_hit = SubElement(search_result, 'search_hit')
            search_hit.set('hit_rank', str(oPepScores.liRanks[iScoreId]))
            search_hit.set('peptide', oPepScores.sOriginalPeptide[1:-1])
            lProteins = oPepScores.sProteinNames.split(',')
            search_hit.set('protein', lProteins[0])
            search_hit.set('num_tot_proteins', str(len(lProteins)))
            search_hit.set('calc_neutral_pep_mass', str(oPepScores.fCalculatedParentMass))
            search_hit.set('massdiff', str(psm.fMeasuredParentMass - oPepScores.fCalculatedParentMass))
            search_hit.set('num_tol_term', '2')
            search_hit.set('num_missed_cleavages', get_num_missed_cleavages(oPepScores.sIdentifiedPeptide, cleave_residues_str))
            # alternative_protein
            if len(lProteins) > 1:
                for iProteins in range(1, len(lProteins)):
                    alternative_protein = SubElement(search_hit, 'alternative_protein')
                    alternative_protein.set('protein', lProteins[iProteins])
            # modification_info
            modification_dict = get_modification_info(oPepScores.sIdentifiedPeptide[1:-1], modification_label_dict)
            if modification_dict:
                modification_info = SubElement(search_hit, "modification_info")
                for key, value in modification_dict.iteritems():
                    mod_aminoacid_mass = SubElement(modification_info, 'mod_aminoacid_mass')
                    mod_aminoacid_mass.set('position', key)
                    mod_aminoacid_mass.set('mass', value)
            # search_score
            search_score = SubElement(search_hit, 'search_score')
            search_score.set('name', score_list_str[iScoreId])
            search_score.set('value', str(oPepScores.lfScores[iScoreId]))

            if iScoreId == 1:
                search_score = SubElement(search_hit, 'search_score')
                search_score.set('name', 'deltacn')
                search_score.set('value', '0.000')
                search_score = SubElement(search_hit, 'search_score')
                search_score.set('name', 'deltacnstar')
                search_score.set('value', '0.000')
                search_score = SubElement(search_hit, 'search_score')
                search_score.set('name', 'spscore')
                search_score.set('value', '0.000')
                search_score = SubElement(search_hit, 'search_score')
                search_score.set('name', 'sprank')
                search_score.set('value', str(oPepScores.liRanks[iScoreId]))
                search_score = SubElement(search_hit, 'search_score')
                search_score.set('name', 'expect')
                search_score.set('value', '0.000')

        iNumProcessedScans += 1
        if iNumProcessedScans % 100 == 0:
            print " Processed & Saved %i Scans\r" % iNumProcessedScans

    if input_filename != '':
        document = ElementTree(root)
        document.write((output_folder + input_filename + '.pep.xml'),
                       encoding='ISO-8859-1',
                       xml_declaration=True,
                       pretty_print=True)

    print 'Done.'
