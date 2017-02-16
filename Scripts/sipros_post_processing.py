'''
Created on Sep 7, 2016

@author: xgo
'''

import getopt, sys, os
import numpy as np
import csv
import math
import re
from sets import Set
from collections import namedtuple
from sklearn import linear_model
from sklearn import preprocessing
# from sklearn.neural_network import MLPClassifier
from subprocess import call
from multiprocessing import Process
from multiprocessing import Queue, cpu_count
import scipy.stats as spst

# # Import Sipros package modules
import sipros_post_module
import sipros_peptides_assembling
import parseconfig

# # Class for ignoring comments '#' in sipros file
CommentedFile = sipros_post_module.CommentedFile

#feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'IPSC', 'OPSC', 'UPSC', 'SPSC', 'pep_psm', 'pro_pep']
feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'IPSC', 'OPSC', 'UPSC', 'SPSC', 'MassWindow', 'PPC', 'OPSC_U', 'OPSC_Ma', 'OPSC_D_Mi', 'PSC_U', 'PSC_Ma', 'PSC_D_Mi']
feature_selection_list = [0, 1, 2, 3, 4, 5]
ptm_str = ['~', '!', '@', '>', '<', '%', '^', '&', '*', '(', ')', '/', '$']
ptm_selection_list = [0]
# ptm_selection_list = [0, 2, 3, 4]
# ptm_selection_list = [0, 1, 2, 3, 4, 5, 6, 7, 8]
for x in ptm_selection_list:
    feature_name_list.append(ptm_str[x])
    

rev_str = 'Rev_' # 'Rev_'
shu_str = 'TestRev_' # 'Shu_'
reserve_str = 'Rev_2_'

rev_str = 'Rev1_' # 'Rev_'
shu_str = 'Rev2_' # 'Shu_'
reserve_str = ''

#rev_str = 'Sh1_' # 'Rev_'
#shu_str = 'Sh2_' # 'Shu_'

#rev_str = 'Rev_'
#shu_str = 'Dec_'

LabelFwd = 1
LabelRev = 2
LabelShu = 3
LabelReserve = 4

mass_tolerance = 0.03

# # Class for PepOutFields object
class PsmFields4(namedtuple('PsmFields',
        ['FileName',  # 0
         'ScanNumber',  # 1
         'ParentCharge',  # 2
         'MeasuredParentMass',  # 3
         'ScanType',  # 4
         'SearchName',  # 5
         'IdentifiedPeptide',  # 6
         'OriginalPeptide',  # 7
         'CalculatedParentMass',  # 8
         'MVH',  # 9
         'Xcorr',  # 10
         'WDP',  # 11
         'ProteinNames',  # 12
         'ScoreAgreement', # 13
         'DeltaRP1',
         'DeltaRP2',
         'DeltaRP3',
         'DeltaRS1',
         'DeltaRS2',
         'DeltaRS3',
         'DiffRP1',
         'DiffRP2',
         'DiffRP3',
         'DiffRS1',
         'DiffRS2',
         'DiffRS3',
         'DiffNorRP1',
         'DiffNorRP2',
         'DiffNorRP3',
         'DiffNorRS1',
         'DiffNorRS2',
         'DiffNorRS3',
         'RetentionTime',
         'Rank',
         'DeltaP'])): # 33 ,
    def __init__(self):
        self.data = self

LabelUnknown = 2
LabelPositive = 1
LabelNegative = 0
LabelFiltered = 3

class PSM:

    iNumScores = 3
    fNeutronMass = 1.00867108694132 # it is Neutron mass
    pattern = re.compile('[^\w\[\]]')

    def __init__(self, psm_field):
        if type(psm_field).__name__ == 'PsmFields3':
            self.pep_psm = int(psm_field.pep_psm)
            self.pro_pep = int(psm_field.pro_pep)
        else:
            self.pep_psm = -1
            self.pro_pep = -1
        self.FileName = psm_field.FileName
        self.bFileNameChanged = False
        self.ScanNumber = int(psm_field.ScanNumber)
        self.ParentCharge = int(psm_field.ParentCharge)
        self.ScanType = psm_field.ScanType
        self.SearchName = psm_field.SearchName
        self.lfScores = [float(psm_field.MVH), float(psm_field.Xcorr), float(psm_field.WDP)]
        self.ProteinNames = psm_field.ProteinNames.strip()
        self.ScoreAgreement = int(psm_field.ScoreAgreement)
        self.IdentifiedPeptide = psm_field.IdentifiedPeptide
        self.OriginalPeptide = psm_field.OriginalPeptide
        self.OriginalPeptide = PSM.pattern.sub('', self.IdentifiedPeptide)
        self.protein_list = []
        self.RealLabel = protein_type(self.ProteinNames, self.protein_list)
        self.lRanks = []
        self.TrainingLabel = LabelNegative  # 0: negative 1: positive 2: unknown
        self.fRankProduct = 0.0
        self.iInnerId = 0
        self.fPredictProbability = 0.0
        self.fMassDiff = 0.0
        self.MeasuredParentMass = float(psm_field.MeasuredParentMass)
        self.CalculatedParentMass = float(psm_field.CalculatedParentMass)
        self.iMassWindow = 0
        self.set_mass_diff(self.MeasuredParentMass, self.CalculatedParentMass)
        
        self.score_differential_list = []
        self.sRTime = '-1.000'
        self.fRtMeasured = 0.0
        self.fRtPredict = 0.0
        self.fRtPvalue = 0.0
        self.iLocalRank = 0
        self.DeltaP = 'NA'
        if type(psm_field).__name__ == 'PsmFields3':
            self.score_differential_list.extend(float(i) for i in psm_field[14:-2])
        elif type(psm_field).__name__ == 'PsmFields4':
            self.score_differential_list.extend(float(i) for i in psm_field[14:32])
            self.sRTime = psm_field.RetentionTime
            self.fRtMeasured = float(self.sRTime)
            self.DeltaP = psm_field.DeltaP 
            self.iLocalRank = int(psm_field.Rank)
        else:
            self.score_differential_list.extend(float(i) for i in psm_field[14:])
        
        if len(self.score_differential_list) != 18:
            print 'error score'
        
        self.NMC = 0
        self.IPSC = 0
        self.OPSC = 0
        self.UPSC = 0 # unique peptide
        self.SPSC = 0 # shared peptide
        self.NRS = 0
        self.PPC = 0
        self.UPPC = 0
        self.SPPC = 0
        self.feature_list = []
        
        self.ML_feature = []
        self.FDR_feature = []
        self.fdr_product = 0
        
        self.OPSC_UMD = [0, 0, 0]
        
        self.SPSC_UMD = [0, 0, 0]
        
        
    def get_feature_list(self):
        del self.feature_list[:]
        '''
        if self.ParentCharge <= 2:
            self.feature_list.append(2)
        else:
            self.feature_list.append(3)
        '''
        self.feature_list.append(self.ParentCharge) # 1: 0
        self.feature_list.extend(self.lfScores) # 2, 3, 4: 1, 2, 3
        self.feature_list.append(self.ScoreAgreement) # 5: 4
        self.feature_list.append(abs(self.fMassDiff)) # 6: 5
        self.feature_list.extend(self.score_differential_list) # 7 - 24: 6 - 23
        self.feature_list.append(self.NMC) # 25: 24
        self.feature_list.append((self.IPSC)) # 26: 25
        # self.OPSC = self.OPSC_UMD[0] + self.OPSC_UMD[1]*0 + self.OPSC_UMD[2]*0 - 1.0
        # self.OPSC = self.OPSC_UMD[0]*0.3333 + self.OPSC_UMD[1]*0.5 + self.OPSC_UMD[2] - 1.0
        self.feature_list.append((self.OPSC)) # 27: 26
        self.feature_list.append((self.UPSC)) # 28: 27
        # self.SPSC = self.SPSC_UMD[0] + self.SPSC_UMD[1]*0 + self.SPSC_UMD[2]*0 - 1.0
        # self.SPSC = self.SPSC_UMD[0]*0.3333 + self.SPSC_UMD[1]*0.5 + self.SPSC_UMD[2] - 1.0
        self.feature_list.append((self.SPSC)) # 29: 28
        '''
        if self.pep_psm != -1:
            self.feature_list.append(self.pep_psm)
            self.feature_list.append(self.pro_pep)
        '''   
        self.feature_list.append(abs(self.iMassWindow)) # 30: 29
        
        # num replicate spectra
        # self.feature_list.append(self.NRS) # 31: 30
        self.feature_list.append((self.PPC)) # 31: 30
        # self.feature_list.append((self.UPPC)) # 31: 30
        # self.feature_list.append((self.SPPC)) # 32: 31
        self.feature_list.extend(self.OPSC_UMD) # 32 - 35: 31 - 34
        self.feature_list.extend(self.SPSC_UMD) # 32 - 35: 31 - 34
        
        for c in ptm_selection_list:
            self.feature_list.append(self.IdentifiedPeptide.count(ptm_str[c])) # 32: 31
        
    
    def set_protein_names(self):
        self.ProteinNames = '{' + ','.join(self.protein_list) + '}'

    def set_feature(self, feature_list):
        del feature_list[:]
        #feature_list.append(self.ParentCharge)
        feature_list.extend(self.lfScores)
        
    def set_mass_diff(self, measured_mass, calculated_mass):
        MassDiffOriginal = measured_mass - calculated_mass
        MassDiff = MassDiffOriginal
        for i in range(-4, 4):
            if abs(MassDiffOriginal - i*PSM.fNeutronMass) < abs(MassDiff):
                MassDiff = MassDiffOriginal - i*PSM.fNeutronMass
                self.iMassWindow = i
        self.fMassDiff = MassDiff
        
    def set_fdr_product(self):
        val = 1.0
        for x in self.FDR_feature:
            val *= (1.0 - x)
        val = 1.0 - pow(val, 1.0/float(len(self.FDR_feature)))
        self.fdr_product = val
        
    def clean_protein_name(self):
        self.ProteinNames = ""
        l = []
        for sProtein in self.protein_list:
            sProtein.strip()
            if not (sProtein.startswith(rev_str)):
                if sProtein not in l:
                    l.append(sProtein)
        self.ProteinNames = '{'+','.join(l) + '}'
        self.protein_list = l

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python xxx.py [options]

Inputs:
    input yyy
    output zzz

Options:
    -h/--help
    -v/--version

Outputs:
    output zzz
'''

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:",
                                    ["help",
                                     "version",
                                     "input",
                                     "config",
                                     "output"])

    # Default working dir and config file
    input_file = ""
    output_folder = ""
    config_file = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print help_message
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print "xxx.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input"):
            input_file = value
        if option in ("-o", "--output"):
            output_folder = value
        if option in ("-c", "--config"):
            config_file = value

    if input_file == "" or output_folder == "":
        print help_message
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')

    return (input_file, config_file, output_folder)

# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    asProteins = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        for sProtein in asProteins:
            sProtein = sProtein.strip()
            if sProtein not in lProtein:
                lProtein.extend(asProteins[:])
    for sProtein in asProteins:
        if not (sProtein.startswith(rev_str) or sProtein.startswith(shu_str)):
            return LabelFwd
    for sProtein in asProteins:
        if sProtein.startswith(shu_str):
            return LabelShu
    if reserve_str != '':
        for sProtein in asProteins:
            if sProtein.startswith(reserve_str):
                return LabelReserve
    return LabelRev

# # split based on score agreement
def categorize_score_agreement(val):
    if val == 3:
        return '3'
    elif val == 2:
        return '2'
    else:
        return '1'

one_top_list = [0, 1, 2, 4]
# # split based on score agreement
def categorize_score_agreement_more_info(val):
    if val in one_top_list:
        return '1'
    elif val == 7:
        return 'A'
    elif val == 3:
        return 'XW'
    elif val == 5:
        return 'MW'
    elif val == 6:
        return 'MX'
    else:
        die("error")

two_top_list = [3, 5, 6]
# # split based on score agreement
def categorize_score_agreement_one_score_agreed(val):
    if val in two_top_list:
        return '2'
    elif val == 7:
        return '3'
    elif val == 0:
        return 'N'
    elif val == 1:
        return 'W'
    elif val == 2:
        return 'X'
    elif val == 4:
        return 'M'
    else:
        die("error")

# # split based on charge
def categorize_parent_charge(val):
    if val == 1:
        return '1'
    elif val == 2:
        return '2'
    elif val >= 3:
        return '3'

# # read the psm table
def read_psm_table(input_file):

    psm_list = []
    
    # read line with csv
    f = csv.reader(CommentedFile(open(input_file, 'rb')),
                            delimiter='\t')
    # skip header
    _sHeader = f.next()
    # get data
    for sLine in f:
        PsmFields_obj = PsmFields4._make(sLine)
        psm_obj = PSM(PsmFields_obj)
        psm_list.append(psm_obj)
        
    i = 0
    for oPsm in psm_list:
        oPsm.iInnerId = i
        i += 1
        
    return (psm_list)

# # Division error handling
divide = sipros_post_module.divide
FDR_parameter = 1.0

# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(TP)
    FDR_accept = True

    if  FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))


# using the global rank
def get_cutoff_rank_product2(FDR_threshold, lPsm):
    F_list = []
    T_list = []
    TT_list = []
    S_list = []
    for oPsm in lPsm:
        if oPsm.RealLabel != LabelRev:
            T_list.append(oPsm.fRankProduct)
            if oPsm.RealLabel == LabelShu:
                S_list.append(oPsm.fRankProduct)
            else:
                TT_list.append(oPsm.fRankProduct)
        else:
            F_list.append(oPsm.fRankProduct)
    F_list = np.array(F_list)
    T_list = np.array(T_list)
    TT_list = np.array(TT_list)
    S_list = np.array(S_list)
    prev_TF_num = 0
    final_cutoff_score = np.amin(T_list, axis=0)
    _final_accept = False
    T_num_best = 0
    TT_num_best = 0
    F_num_best = 0
    S_num_best = 0
    TF_num_best = 0
    for cutoff_score in T_list:
        F_num = (F_list <= cutoff_score).sum()
        T_num = (T_list <= cutoff_score).sum()
        TT_num = (TT_list <= cutoff_score).sum()
        S_num = (S_list <= cutoff_score).sum()
        # (F_num, T_num) = get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex)
        TF_num = F_num + T_num
        TF_num = S_num + TT_num
        (FDR_accept, FDR_value) = FDR_calculator(S_num, TT_num)
        # update final values if conditions satisfies
        # 1) FDR_accept is True
        # 2) FDR_value should be less than or equal to FDR_threshold_1
        # 3) TF_num is greater than to previous TF_num
        if (FDR_accept is True) and (FDR_value <= FDR_threshold) and (TF_num > prev_TF_num) :
            final_cutoff_score = cutoff_score
            _final_accept = FDR_accept
            # previous TF_num
            prev_TF_num = TF_num
            T_num_best = T_num
            F_num_best = F_num
            S_num_best = S_num
            TF_num_best = TF_num

    print "Q-value_cutoff:\t%f\t%d\t%d\t%d" % (final_cutoff_score, T_num_best - S_num_best, F_num_best, S_num_best)
    return final_cutoff_score


# using the q value rank
def get_cutoff_q_rank_product(FDR_threshold, lPsm):
    iNumScores = lPsm[0].iNumScores

    T_num = 0
    F_num = 0
    S_num = 0
    for oPsm in lPsm:
        if oPsm.RealLabel != LabelRev:
            T_num += 1
            if oPsm.RealLabel == LabelShu:
                S_num += 1
        else:
            F_num += 1
    print "Before Filtering:\t\t%d\t%d\t%d" % (T_num - S_num, F_num, S_num)

    for i in range(iNumScores):
        newlist = sorted(lPsm, key=lambda x: x.lfScores[i], reverse=True)
        T_num = 0
        F_num = 0
        for j in range(len(newlist)):
            if newlist[j].RealLabel != LabelRev:
                T_num += 1
            else:
                F_num += 1
            (_FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)
            newlist[j].lRanks.append(FDR_value)
        fSmallestQ = 1
        for j in range(len(newlist) - 1, -1, -1):
            if fSmallestQ > newlist[j].lRanks[i]:
                fSmallestQ = newlist[j].lRanks[i]
            if newlist[j].lRanks[i] > fSmallestQ:
                newlist[j].lRanks[i] = fSmallestQ

    for oPsm in lPsm:
        # fTemp = oPsm.lRanks[0] * oPsm.lRanks[1] * oPsm.lRanks[2]
        # oPsm.fRankProduct = np.power(float(fTemp), 1.0 / 3.0)
        fTemp = (1 - ((1 - oPsm.lRanks[0]) * (1 - oPsm.lRanks[1]) * (1 - oPsm.lRanks[2])))
        oPsm.fRankProduct = fTemp

    final_cutoff_score = get_cutoff_rank_product2(FDR_threshold, lPsm)
    return final_cutoff_score

def show_TP_TN_FP_FN(label_np, predict_np):
    true_np = (label_np == predict_np)
    TP = label_np[true_np].sum()
    TN = (label_np[true_np] == 0).sum()
    false_np = (label_np != predict_np)
    FP = (label_np[false_np] == 0).sum()
    FN = label_np[false_np].sum()
    print "TP\tTN\tFP\tFN"
    print "%d\t%d\t%d\t%d" % (TP, TN, FP, FN)
    
def show_Fdr_category(psm_dict):
    Best_last_list = [0, 0, 0]
    for _key, lPsm in psm_dict.iteritems():
        #print _key
        list_sorted = sorted(lPsm, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
        T_num = 0
        F_num = 0
        Fwd_num = 0
        Rev_num = 0
        Shu_num = 0
        Best_list = [0, 0, 0]
        for oPsm in list_sorted:
            if oPsm.RealLabel == LabelFwd:
                T_num += 1
                Fwd_num += 1
            elif oPsm.RealLabel == LabelRev:
                F_num += 1
                Rev_num += 1
            else:
                T_num += 1
                Shu_num += 1
            (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
            if (FDR_accept is True) and (FDR_value <= 0.01): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
                if Best_list[0] < Fwd_num:
                    Best_list = [Fwd_num, Rev_num, Shu_num]
        #sys.stdout.write('%d\t%d\t%d\n' % (Best_list[0], Best_list[1], Best_list[2]))
        Best_last_list = [Best_last_list[0] + Best_list[0], Best_last_list[1] + Best_list[1], Best_last_list[2] + Best_list[2]]
    sys.stdout.write('%d\t%d\t%d\t' % (Best_last_list[0], Best_last_list[1], Best_last_list[2]))
    pass

def filter_Fdr(psm_list, fdr):
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    psm_filtered_list = []
    psm_left_list = []
    cutoff_probability = 0.0
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
    
    sys.stdout.write('\n%d\t%d\t%d\t\n' % (Best_list[0], Best_list[1], Best_list[2]))
    
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
            oPsm.TrainingLabel = LabelFiltered
        else:
            psm_left_list.append(oPsm)
    
    return (psm_filtered_list, psm_left_list)

def filter_Fdr2(psm_list, fdr):
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    psm_filtered_list = []
    psm_left_list = []
    cutoff_probability = 0.0
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fRtPvalue < 0.0001:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
    
    sys.stdout.write('\n%d\t%d\t%d\t\n' % (Best_list[0], Best_list[1], Best_list[2]))
    
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
            oPsm.TrainingLabel = LabelFiltered
        else:
            psm_left_list.append(oPsm)
    
    return (psm_filtered_list, psm_left_list)

def show_Fdr_charge(psm_list):
    psm_new_list_1 = []
    psm_new_list_2 = []
    for oPsm in psm_list:
        if oPsm.ParentCharge <= 2:
            psm_new_list_1.append(oPsm)
        else:
            psm_new_list_2.append(oPsm)
    show_Fdr(psm_new_list_1, None, None)
    show_Fdr(psm_new_list_2, None, None)

psm_num_list = []
pep_num_list = []

def show_Fdr_varied(psm_list, fdr):
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]

    psm_filtered_list = []
    
    cutoff_probability = 0.0
    
    peptide_set = Set()
    num_fwr_pep = 0
    num_shu_pep = 0
    best_fwr_pep = 0
    best_shu_pep = 0
    
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelShu:
                num_shu_pep += 1
                peptide_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        elif oPsm.RealLabel == LabelShu:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
        (FDR_accept, FDR_value) = FDR_calculator(num_shu_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (best_fwr_pep + best_shu_pep) < (num_fwr_pep + num_shu_pep):        
                best_fwr_pep = num_fwr_pep
                best_shu_pep = num_shu_pep
            
    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
    
    # sys.stdout.write('\n'+str(len(psm_filtered_list))+'\t')
    
    print "{:,d} ({:.2f}%)\t{:,d} ({:.2f}%)".format(Best_list[0], (100* float(Best_list[2])/float(Best_list[0])), best_fwr_pep, (100.0*float(best_shu_pep)/float(best_fwr_pep)))
    psm_num_list.append(Best_list[0])
    pep_num_list.append(best_fwr_pep)
    return psm_filtered_list

def show_Fdr(psm_list, sKey, fdr=None):
    
    fdr_float = 0.01
    if fdr !=None:
        fdr_float = fdr
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    '''
    fwr = open("fwr.txt", "a")
    rev = open("rev.txt", "a")
    '''
    
    peptide_set = Set()
    num_fwr_pep = 0
    num_shu_pep = 0
    best_fwr_pep = 0
    best_shu_pep = 0
    
    psm_filtered_list = []
    cutoff_probability = 0.0
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelShu:
                num_shu_pep += 1
                peptide_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            Rev_num += 1
        elif oPsm.RealLabel == LabelShu:
            Shu_num += 1
        else:
            print 'error 768.'
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr_float): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
                best_fwr_pep = num_fwr_pep
                best_shu_pep = num_shu_pep
            

    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    sys.stdout.write("{:,d}\n{:.2f}%\n{:.2f}%\n{:,d}\n{:.2f}%\n".format(Best_list[0], 100.0*float(Best_list[1])/float(Best_list[0]), 100.0*float(Best_list[2])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    # sys.stdout.write("%d\t[%.2f%%]\t[%.2f%%]\t%d\t[%.2f%%]\n" % (Best_list[0], 100.0*float(Best_list[2])/float(Best_list[0]), 100.0*float(Best_list[1])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    
    return psm_filtered_list

def show_Fdr_Pep(psm_list, sKey, fdr=None):
    
    fdr_float = 0.01
    if fdr !=None:
        fdr_float = fdr
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    '''
    fwr = open("fwr.txt", "a")
    rev = open("rev.txt", "a")
    '''
    
    peptide_set = Set()
    num_fwr_pep = 0
    num_shu_pep = 0
    best_fwr_pep = 0
    best_shu_pep = 0
    
    psm_filtered_list = []
    cutoff_probability = 0.0
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelShu:
                num_shu_pep += 1
                peptide_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        elif oPsm.RealLabel == LabelShu:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(num_shu_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr_float): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
                best_fwr_pep = num_fwr_pep
                best_shu_pep = num_shu_pep
            

    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    
    sys.stdout.write("%d\t[%.2f%%]\t[%.2f%%]\t%d\t[%.2f%%]\n" % (Best_list[0], 100.0*float(Best_list[2])/float(Best_list[0]), 100.0*float(Best_list[1])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    
    return psm_filtered_list

def show_Fdr_peptide(psm_list, sKey, fdr=None):
    
    fdr_float = 0.01
    if fdr !=None:
        fdr_float = fdr
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]

    psm_filtered_list = []
    cutoff_probability = 0.0
    peptide_set = Set()
    # without considering training label
    for oPsm in list_sorted:
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str in peptide_set:
            continue
        else:
            peptide_set.add(pep_str)
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr_float): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
            

    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
   
    sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    
    sys.stdout.write('%d\t%d\t%d\t' % (Best_list[0], Best_list[1], Best_list[2]))
    
    return psm_filtered_list

def logistic_regression(psm_dict, psm_list):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    psm_list_selected = []
    for key, lPsm in psm_dict.iteritems():
        sys.stdout.write(key + "\t")
        psm_list_selected = lPsm
        data_list = []
        label_list = []
        unknown_list = []
        id_list = []
        for oPsm in lPsm:
            #feature_list = []
            #oPsm.set_feature(feature_list)
            if len(oPsm.feature_list) == 0:
                print 'check'
            unknown_list.append(oPsm.feature_list)
            id_list.append(oPsm.iInnerId)
            if oPsm.TrainingLabel != LabelUnknown:
                data_list.append(oPsm.feature_list)
                label_list.append(oPsm.TrainingLabel)

        data_np = np.array(data_list)
        sys.stdout.write(str(len(data_list)) + "\t")
        label_np = np.array(label_list)
        
        '''
        for i in range(len(feature_name_list)):
            del feature_selection_list[:]
            feature_selection_list.append(i)
            for j in range(len(feature_name_list)):
                if j != i:
                    feature_selection_list.append(j)
        '''    
        '''
        del feature_selection_list[:]
        for i in range(len(feature_name_list)):
            if i < 6 or i > 23 or ( i>=15 and i<=17):
                feature_selection_list.append(i)
        #feature_selection_list.extend([1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])
        '''
        train_data_np = data_np[:, feature_selection_list]
    
    # # training
        class_dict = {0: 1, 1:1000}
        logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
        logreg.fit(train_data_np, label_np)
        predict_np = logreg.predict(train_data_np)
        #np.savetxt("xxx.txt", train_data_np)
        #return

    #show_TP_TN_FP_FN(label_np, predict_np)
    #print 'Actual number of iterations for all classes.'
        #print logreg.n_iter_
    # # test
        unknown_np = np.array(unknown_list)
        test_unknown_np = unknown_np[:, feature_selection_list]
        predict_np = logreg.predict_proba(test_unknown_np)

        for i in range(len(predict_np)):
            psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]

    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
        psm_filtered_list_local = show_Fdr(psm_list_selected, key)
        psm_filtered_list.extend(psm_filtered_list_local)
    
    #print 'Coefficient of the features in the decision function:'
    #print logreg.coef_
        idx = 0
        for i in range(len(feature_name_list)):
            if i in feature_selection_list:
                sys.stdout.write('%.3f' % logreg.coef_[0][idx])
                idx += 1
            sys.stdout.write('\t')
        sys.stdout.write('\n')
        #print str(logreg.intercept_)
    '''
    for x in logreg.coef_[0]:
        sys.stdout.write(str(x))
        sys.stdout.write('\t')
    '''
    return psm_filtered_list

def re_rank(psm_list):
    psm_new_list = []
    psm_dict = {}
    for oPsm in psm_list:
        sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
        if sId in psm_dict:
            if oPsm.fPredictProbability > psm_dict[sId].fPredictProbability:
                psm_dict[sId] = oPsm
        else:
            psm_dict[sId] = oPsm
    
    for _key, value in psm_dict.iteritems():
        psm_new_list.append(value)
    
    return psm_new_list

def writeout_feature(list_before_filtering, list_after_filtering):
    filename_before_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/fig2/fig2_Angelo_before.txt'
    filename_after_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/fig2/fig2_Angelo_after.txt'
    ms2_filename_dict = {}
    ms2_filename_str = ''
    index_int = 0
    scan_set = Set()
    id_str = ''
    with open(filename_before_str, 'w') as fw:
        fw.write('RealLabel\tMVH\tXcorr\tWDP\tDiffMVH\tDiffXcorr\tDiffWDP\tMassError\tMissedCleavage\tOPSC\tSPSC\n')
        for oPsm in list_before_filtering:
            if oPsm.FileName in ms2_filename_dict:
                ms2_filename_str = ms2_filename_dict[oPsm.FileName]
            else:
                ms2_filename_dict[oPsm.FileName] = str(index_int)
                ms2_filename_str = str(index_int)
                index_int += 1
            id_str = ms2_filename_str + '_' + str(oPsm.ScanNumber) + '_' + str(oPsm.ParentCharge) + '_' + oPsm.IdentifiedPeptide
            if id_str in scan_set:
                continue
            else:
                scan_set.add(id_str)
                fw.write(str(oPsm.RealLabel))
                fw.write('\t')
                fw.write(str(oPsm.lfScores[0]))
                fw.write('\t')
                fw.write(str(oPsm.lfScores[1]))
                fw.write('\t')
                fw.write(str(oPsm.lfScores[2]))
                fw.write('\t')
                fw.write(str(oPsm.score_differential_list[9]))
                fw.write('\t')
                fw.write(str(oPsm.score_differential_list[10]))
                fw.write('\t')
                fw.write(str(oPsm.score_differential_list[11]))
                # fw.write('\t')
                # fw.write(str(abs(oPsm.iMassWindow)))
                fw.write('\t')
                fw.write(str(abs(oPsm.fMassDiff)))
                fw.write('\t')
                fw.write(str(oPsm.NMC))
                # fw.write('\t')
                # fw.write(str(oPsm.IPSC))
                fw.write('\t')
                fw.write(str(oPsm.OPSC))
                # fw.write('\t')
                # fw.write(str(oPsm.UPSC))
                fw.write('\t')
                fw.write(str(oPsm.SPSC))
                fw.write('\n')
    
    with open(filename_after_str, 'w') as fw:
        fw.write('RealLabel\tMVH\tXcorr\tWDP\tDiffMVH\tDiffXcorr\tDiffWDP\tMassError\tMissedCleavage\tOPSC\tSPSC\n')
        for oPsm in list_after_filtering:
            if oPsm.FileName in ms2_filename_dict:
                ms2_filename_str = ms2_filename_dict[oPsm.FileName]
            else:
                ms2_filename_dict[oPsm.FileName] = str(index_int)
                ms2_filename_str = str(index_int)
                index_int += 1
            fw.write(str(oPsm.RealLabel))
            fw.write('\t')
            fw.write(str(oPsm.lfScores[0]))
            fw.write('\t')
            fw.write(str(oPsm.lfScores[1]))
            fw.write('\t')
            fw.write(str(oPsm.lfScores[2]))
            fw.write('\t')
            fw.write(str(oPsm.score_differential_list[9]))
            fw.write('\t')
            fw.write(str(oPsm.score_differential_list[10]))
            fw.write('\t')
            fw.write(str(oPsm.score_differential_list[11]))
            fw.write('\t')
            fw.write(str(abs(oPsm.iMassWindow)))
            fw.write('\t')
            fw.write(str(abs(oPsm.fMassDiff)))
            fw.write('\t')
            fw.write(str(oPsm.NMC))
            # fw.write('\t')
            # fw.write(str(oPsm.IPSC))
            fw.write('\t')
            fw.write(str(oPsm.OPSC))
            # fw.write('\t')
            # fw.write(str(oPsm.UPSC))
            fw.write('\t')
            fw.write(str(oPsm.SPSC))
            fw.write('\n')

def logistic_regression_no_category(psm_dict, psm_list, psm_neg_list, fdr_given=None, output_folder=None, input_file=None, config_dict=None):
    # machine learning
    # # construct training data
    #psm_list_selected = []
    train_data_list = []
    train_label_list = []
    test_data_list = []
    '''
    for key, lPsm in psm_dict.iteritems():
        #sys.stdout.write(key + "\t")
        for oPsm in lPsm:
            if len(oPsm.feature_list) == 0:
                print 'check'
            if oPsm.TrainingLabel == LabelFiltered:
                continue
            test_data_list.append(oPsm.feature_list)
            if oPsm.TrainingLabel != LabelUnknown:
                train_data_list.append(oPsm.feature_list)
                train_label_list.append(oPsm.TrainingLabel)
    '''
    psm_rank_list = []
    '''
    psm_rank_U_list = []
    psm_rank_M_list = []
    psm_rank_D_list = []
    '''
    num_feature_int = (31 + len(ptm_selection_list))
    positive_int = 1
    negative_int = 0
    bDisableLocalRank = False
    if len(psm_list) < 800000:
        bDisableLocalRank = True
    for oPsm in psm_list:
        if len(oPsm.feature_list) != num_feature_int:
            pass
            # print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelReserve:
            # train_data_list.append(oPsm.feature_list)
            # train_label_list.append(negative_int)
            continue
        # if oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1:
        test_data_list.append(oPsm.feature_list)
        psm_rank_list.append(oPsm)
        '''
        if oPsm.iLocalRank == 0:
            psm_rank_U_list.append(oPsm)
        elif oPsm.iLocalRank == 2:
            psm_rank_D_list.append(oPsm)
        else:
            psm_rank_M_list.append(oPsm)
        '''
        # if oPsm.TrainingLabel != LabelUnknown and oPsm.iLocalRank == 0:
        # if oPsm.RealLabel == LabelFwd and (oPsm.iLocalRank == 0 ):
        if oPsm.RealLabel != LabelRev and (oPsm.iLocalRank == 0 or bDisableLocalRank ):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(positive_int)
        elif oPsm.RealLabel == LabelRev: # and (oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(negative_int)
                
    sys.stdout.write(str(len(train_data_list)) + "\t")
        
    train_data_np = np.array(train_data_list)[:, feature_selection_list]
    train_label_np = np.array(train_label_list)
    
     # only forward left
    unique_np = np.unique(train_label_np)
    if unique_np.shape[0] == 1:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
        return psm_filtered_list_local
    
    # # training
    # num_positive = float((train_label_np==LabelPositive).sum())
    # num_negative = float((train_label_np==LabelNegative).sum())
    # num_positive = 100.0
    # num_negative = 1.0
    # class_weight_dict = {0: (num_positive/(num_negative+num_positive)), 1:(num_negative/(num_negative+num_positive))}
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)

    # # test
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i, 1]
        psm_rank_list[i].ML_feature.append(predict_np[i, 1])
    # fdr_rank(psm_list, 0)
    # return None
    '''
    # U
    test_data_list = []
    for oPsm in psm_rank_U_list:
        test_data_list.append(oPsm.feature_list)
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    for i in range(len(predict_np)):
        psm_rank_U_list[i].fPredictProbability = predict_np[i, 1]
    psm_rank_U_list = re_rank(psm_rank_U_list)
    psm_rank_U_list = show_Fdr(psm_rank_U_list, None, fdr=fdr_given)
    # M
    test_data_list = []
    for oPsm in psm_rank_M_list:
        test_data_list.append(oPsm.feature_list)
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    for i in range(len(predict_np)):
        psm_rank_M_list[i].fPredictProbability = predict_np[i, 1]
    psm_rank_M_list = re_rank(psm_rank_M_list)
    psm_rank_M_list = show_Fdr(psm_rank_M_list, None, fdr=fdr_given)
    # D
    test_data_list = []
    for oPsm in psm_rank_D_list:
        test_data_list.append(oPsm.feature_list)
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    for i in range(len(predict_np)):
        psm_rank_D_list[i].fPredictProbability = predict_np[i, 1]
    psm_rank_D_list = re_rank(psm_rank_D_list)
    psm_rank_D_list = show_Fdr(psm_rank_D_list, None, fdr=fdr_given)
    psm_rank_list = psm_rank_U_list
    psm_rank_list.extend(psm_rank_M_list)
    psm_rank_list.extend(psm_rank_D_list)
    '''
    psm_new_list = re_rank(psm_rank_list)
    # del psm_list[:]
    # psm_list.extend(psm_new_list)
    
    # psm_new_list = mass_filter(psm_new_list)
    
    for fdr_f in [0.005, 0.01, 0.02]:
        temp_list = show_Fdr_varied(psm_new_list, fdr_f/3.0)
        '''
        folder_str = output_folder + 'fdr_' + str(fdr_f) +'/'
        if not os.path.exists(folder_str):
            os.makedirs(folder_str)
        else:
            for the_file in os.listdir(folder_str):
                file_path = os.path.join(folder_str, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                        # elif os.path.isdir(file_path): shutil.rmtree(file_path)
                except Exception as e:
                    print(e)
        generate_psm_pep_txt(input_file, folder_str, temp_list)
        '''
    for num in psm_num_list:
        sys.stdout.write('{:,d}\t'.format(num))
    sys.stdout.write('\n')
    for num in pep_num_list:
        sys.stdout.write('{:,d}\t'.format(num))
    sys.stdout.write('\n')
    
    
    if config_dict != None:
        if config_dict[FDR_Filtering_str] == 'PSM':
            psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
        else:
            psm_filtered_list_local = show_Fdr_Pep(psm_new_list, None, fdr=fdr_given)
        # writeout_feature(psm_list, psm_filtered_list_local)
        return psm_filtered_list_local
    
    psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
        # psm_filtered_list_local = show_Fdr_peptide(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    
    
    
    psm_filtered_list = []
    psm_filtered_list.extend(psm_filtered_list_local)
    b = logreg.coef_[0]
    sum_float = 0
    for val in b:
        sum_float += abs(val)
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            # sys.stdout.write('%.3f' % logreg.coef_[0][idx])
            sys.stdout.write('%.5f' % (b[idx]/sum_float))
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
    '''
    order_list = [2, 0, 1, 6, 4, 5, 12, 3, 7, 8, 9, 10, 11]
    print "\nfeature importance:"
    for ind in order_list:
        sys.stdout.write('%.5f\n' % logreg.coef_[0][ind])
    '''    
        
    '''
    iC = 0
    for oPsm in psm_filtered_list:
        if oPsm.ScoreAgreement == 1 and oPsm.RealLabel == LabelFwd:
            iC += 1
    print str(iC)
    '''        
        
    return psm_filtered_list

def logistic_regression_1_LR(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None, output_folder=None, input_file=None):
    # machine learning
    # # construct training data
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    
    positive_int = 1
    negative_int = 0
    
    test_list_1 = []
    train_list_1 = []
    label_list_1 = []
    psm_list_1 = []
    
    test_list_2 = []
    train_list_2 = []
    label_list_2 = []
    psm_list_2 = []
    
    psm_list_3 = []
    test_list_3 = []
    
    num_feature_int = (31 + len(ptm_selection_list))           
    for oPsm in psm_list:
        if len(oPsm.feature_list) != num_feature_int:
            print 'check'
        if oPsm.iLocalRank == 0: # U
            test_list_1.append(oPsm.feature_list)
            psm_list_1.append(oPsm)
            train_list_1.append(oPsm.feature_list)
            if oPsm.RealLabel != LabelRev:
                label_list_1.append(positive_int)
                # train_list_2.append(oPsm.feature_list)
                # label_list_2.append(positive_int)
            else:
                label_list_1.append(negative_int)
        else:
            if oPsm.iLocalRank == 1:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRev:
                    # pass
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(positive_int)
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            elif oPsm.iLocalRank == 3:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRev:
                    pass
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            else:
                test_list_3.append(oPsm.feature_list)
                psm_list_3.append(oPsm)
                if oPsm.RealLabel != LabelRev:
                    pass
                else:
                    pass
                    # train_list_1.append(oPsm.feature_list)
                    # label_list_1.append(negative_int)
    
    # # training
    print 'train 1 size: %d' % len(train_list_1)
    train_data_np = np.array(train_list_1)[:, feature_selection_list]
    train_label_np = np.array(label_list_1)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_1)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_1[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_1)
    psm_list_1 = show_Fdr(psm_new_list, None, fdr=fdr_given)
    
    test_data_np = np.array(test_list_3)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_3[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_3)
    psm_list_3 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
        
    # # training
    print 'train 2 size: %d' % len(train_list_2)
    train_data_np = np.array(train_list_2)[:, feature_selection_list]
    train_label_np = np.array(label_list_2)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_2)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_2[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_2)
    psm_list_2 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
    
    psm_list_all = []
    psm_list_all.extend(psm_list_1)
    psm_list_all.extend(psm_list_2)
    psm_list_all.extend(psm_list_3)
    psm_new_list = re_rank(psm_list_all)
    psm_list_all = show_Fdr(psm_list_all, None, fdr=fdr_given)


def logistic_regression_2_LR(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None, output_folder=None, input_file=None):
    # machine learning
    # # construct training data
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    
    positive_int = 1
    negative_int = 0
    
    test_list_1 = []
    train_list_1 = []
    label_list_1 = []
    psm_list_1 = []
    
    test_list_2 = []
    train_list_2 = []
    label_list_2 = []
    psm_list_2 = []
    
    psm_list_3 = []
    test_list_3 = []
    
    num_feature_int = (31 + len(ptm_selection_list))           
    for oPsm in psm_list:
        if len(oPsm.feature_list) != num_feature_int:
            # pass
            print 'check'
        if oPsm.iLocalRank == 0: # U
            test_list_1.append(oPsm.feature_list)
            psm_list_1.append(oPsm)
            train_list_1.append(oPsm.feature_list)
            if oPsm.RealLabel != LabelRev:
                label_list_1.append(positive_int)
                # train_list_2.append(oPsm.feature_list)
                # label_list_2.append(positive_int)
            else:
                label_list_1.append(negative_int)
        else:
            if oPsm.iLocalRank == 1:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRev:
                    # pass
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(positive_int)
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            elif oPsm.iLocalRank == 3:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRev:
                    pass
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            else:
                test_list_3.append(oPsm.feature_list)
                psm_list_3.append(oPsm)
                if oPsm.RealLabel != LabelRev:
                    pass
                else:
                    pass
                    # train_list_1.append(oPsm.feature_list)
                    # label_list_1.append(negative_int)
    
    # # training
    print 'train 1 size: %d' % len(train_list_1)
    train_data_np = np.array(train_list_1)[:, feature_selection_list]
    train_label_np = np.array(label_list_1)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_1)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_1[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_1)
    psm_list_1 = show_Fdr(psm_new_list, None, fdr=fdr_given)
    
    test_data_np = np.array(test_list_3)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_3[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_3)
    psm_list_3 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
        
    # # training
    print 'train 2 size: %d' % len(train_list_2)
    train_data_np = np.array(train_list_2)[:, feature_selection_list]
    train_label_np = np.array(label_list_2)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_2)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_2[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_2)
    psm_list_2 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
    
    psm_list_all = []
    psm_list_all.extend(psm_list_1)
    psm_list_all.extend(psm_list_2)
    psm_list_all.extend(psm_list_3)
    psm_new_list = re_rank(psm_list_all)
    psm_list_all = show_Fdr(psm_list_all, None, fdr=fdr_given)

import random

def test_train_bias(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    train_data_list = []
    train_label_list = []
    train_id_list = []
    test_data_list = []
    test_id_list = []
    # # psm with the same name
    psm_sa_1_dict = {}
    
    for oPsm in psm_list:
        if oPsm.ScoreAgreement <= 1:
            sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
            if sId in psm_sa_1_dict:
                if psm_sa_1_dict[sId] == True:
                    train_data_list.append(oPsm.feature_list)
                    if oPsm.RealLabel == LabelRev:
                        train_label_list.append(0)
                    else:
                        train_label_list.append(1)
                    train_id_list.append(oPsm.iInnerId)
                else:
                    test_data_list.append(oPsm.feature_list)
                    test_id_list.append(oPsm.iInnerId)
            else:
                f_rand = random.random()
                if f_rand >= 0.5:
                    psm_sa_1_dict[sId] = True
                    train_data_list.append(oPsm.feature_list)
                    if oPsm.RealLabel == LabelRev:
                        train_label_list.append(0)
                    else:
                        train_label_list.append(1)
                    train_id_list.append(oPsm.iInnerId)
                else:
                    psm_sa_1_dict[sId] = False
                    test_data_list.append(oPsm.feature_list)
                    test_id_list.append(oPsm.iInnerId)
        else:
            f_rand = random.random()
            if f_rand >= 0.5:
                train_data_list.append(oPsm.feature_list)
                if oPsm.RealLabel == LabelRev:
                    train_label_list.append(0)
                else:
                    train_label_list.append(1)
                train_id_list.append(oPsm.iInnerId)
            else:
                test_data_list.append(oPsm.feature_list)
                test_id_list.append(oPsm.iInnerId)
                
    train_data_np = np.array(train_data_list)
    sys.stdout.write(str(len(train_data_list)) + "\t")
    train_label_np = np.array(train_label_list)
        
    train_data_np = train_data_np[:, feature_selection_list]
   
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)

    # # test
    test_data_np = np.array(test_data_list)
    test_unknown_np = test_data_np[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[test_id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    print 'testing results:'
    for i in range(1, 11):
        fdr_f = 0.001 * i
        show_Fdr_varied(psm_new_list, fdr_f)
    
    for i in range(len(predict_np)):
        psm_list[test_id_list[i]].fPredictProbability = 0 
    predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[train_id_list[i]].fPredictProbability = predict_np[i, 1]
    print 'train results:'
    psm_new_list = re_rank(psm_list)
    for i in range(1, 11):
        fdr_f = 0.001 * i
        show_Fdr_varied(psm_new_list, fdr_f)       
        
    return None

from sklearn import svm

def svm_one_class_test(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel == LabelNegative:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
    
    clf.fit(train_data_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
        
    return psm_filtered_list

def logistic_regression_no_category_rt(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    for key, lPsm in psm_dict.iteritems():
        #sys.stdout.write(key + "\t")
        for oPsm in lPsm:
            if len(oPsm.feature_list) == 0:
                print 'check'
                
            if oPsm.RealLabel == LabelRev:
                data_list.append(oPsm.feature_list)
                label_list.append(oPsm.TrainingLabel)
            if oPsm.TrainingLabel == LabelFiltered:
                if oPsm.RealLabel != LabelRev:
                    data_list.append(oPsm.feature_list)
                    label_list.append(oPsm.TrainingLabel)
                continue
            unknown_list.append(oPsm.feature_list)
            id_list.append(oPsm.iInnerId)
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
        
    train_data_np = data_np[:, feature_selection_list]
    '''
    if not psm_left_list is None:
        np.savetxt("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/train_data_np.txt", train_data_np)
    '''
    # # training
    #class_dict = {0: 1, 1:1000}
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]

    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, key, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, key, fdr=fdr_given)
    psm_filtered_list.extend(psm_filtered_list_local)
    # sys.stdout.write('\n')
    # show_Fdr_category(psm_dict)
    #print 'Coefficient of the features in the decision function:'
    #print logreg.coef_
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % logreg.coef_[0][idx])
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
        
    return psm_filtered_list

def get_num_missed_cleavage_sites(sIdentifiedSeq, sResiduesBeforeCleavage, sResiduesAfterCleavage):
    count_int = 0
    for i in range(len(sIdentifiedSeq) - 1):
        if sIdentifiedSeq[i] in sResiduesBeforeCleavage and sIdentifiedSeq[i + 1] in sResiduesAfterCleavage:
            count_int += 1
    return count_int

def generate_pep_pro_features(lPsm, config_dict):
    # number of replicate spectra
    identified_peptide_dict = {}
    # number of sibling ions
    identified_peptide_with_different_charges_dict = {}
    # peptide to protein dictionary
    original_pep_pro_dict = {}
    for oPsm in lPsm:
        # number of missed cleavage sites
        '''
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[cleave_after_residues_str], 
                                                 config_dict[cleave_before_residues_str])
        '''
        if oPsm.IdentifiedPeptide in identified_peptide_with_different_charges_dict:
            if oPsm.ParentCharge not in identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide]:
                identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide].append(oPsm.ParentCharge)
        else:
            identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide] = [oPsm.ParentCharge]
        if oPsm.OriginalPeptide not in original_pep_pro_dict:
            original_pep_pro_dict[oPsm.OriginalPeptide] = Set()
            for pro_str in oPsm.protein_list:
                original_pep_pro_dict[oPsm.OriginalPeptide].add(pro_str)
        else:
            for pro_str in oPsm.protein_list:
                original_pep_pro_dict[oPsm.OriginalPeptide].add(pro_str)
        sId = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if sId in identified_peptide_dict:
            identified_peptide_dict[sId] += 1
        else:
            identified_peptide_dict[sId] = 1
        
    # number of sibling modification
    original_peptide_dict = {}
    pattern = re.compile('[\W_]+')
    for (pep_str, charge_list) in identified_peptide_with_different_charges_dict.iteritems():
        # get the original peptide
        pep_original_str = '[' + pattern.sub('', pep_str) + ']'
        if pep_original_str not in original_peptide_dict:
            original_peptide_dict[pep_original_str] = 1
        else:
            original_peptide_dict[pep_original_str] += 1
    # number of sibling ions
    for oPsm in lPsm:
        if oPsm.IdentifiedPeptide not in identified_peptide_with_different_charges_dict:
            print 'error in number of sibling ions'
        else:
            pass
            # oPsm.IPSC = len(identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide])
        # # get the number of sibling modification
        if oPsm.OriginalPeptide not in original_peptide_dict:
            print 'error in number of sibling modification'
        else:
            pass
            # oPsm.OPSC = original_peptide_dict[oPsm.OriginalPeptide]
        # # get the number of replicate spectra
        sId = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if sId in identified_peptide_dict:
            oPsm.NRS = identified_peptide_dict[sId]
        else:
            print 'error in number of replicate spectra'
    
    # protein to peptide dictionary
    pro_original_pep_dict = {}
    for (pep_str, pro_list) in original_pep_pro_dict.iteritems():
        for pro_str in pro_list:
            if pro_str not in pro_original_pep_dict:
                pro_original_pep_dict[pro_str] = Set()
                pro_original_pep_dict[pro_str].add(pep_str)
            else:
                pro_original_pep_dict[pro_str].add(pep_str)
    
    # number of sibling unique peptide and number of sibling shared peptide
    '''
    pro_unique_dict = {}
    pro_shared_dict = {}
    changed_flag = False
    for oPsm in lPsm:
        pro_list = original_pep_pro_dict[oPsm.OriginalPeptide]
        changed_flag = False
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                oPsm.protein_list.append(protein)
                changed_flag = True
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
        if len(pro_list) > 1:
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
        else:
            pro_str = pro_list.pop()
            pro_list.add(pro_str)
            if pro_str in pro_unique_dict:
                pro_unique_dict[pro_str] += 1
            else:
                pro_unique_dict[pro_str] = 1
    
    # collect features
    for oPsm in lPsm:
        pro_list = original_pep_pro_dict[oPsm.OriginalPeptide]
        for protein in pro_list:
            if protein in pro_unique_dict:
                oPsm.UPSC += pro_unique_dict[protein]
                if len(pro_list) == 1:
                    oPsm.UPSC -= 1
            if protein in pro_shared_dict:
                oPsm.SPSC += pro_shared_dict[protein]
                if len(pro_list) > 1:
                    oPsm.SPSC -= 1
    
    '''
    
    # prophet way
    num_unique_per_pro = 0
    num_shared_per_pro = 0
    num_balanced_shared_per_pro = 0
    
    max_unique_per_psm = 0
    max_shared_per_psm = 0
    max_balanced_shared_per_psm = 0
    
    num_unique_per_psm = 0
    num_shared_per_psm = 0
    num_balanced_shared_per_psm = 0
    
    num_per_pro = 0 # linked to a protein
    num_balanced_per_pro = 0 # linked to a protein
    
    max_per_psm = 0 # linked to protein
    max_balanced_per_psm = 0 # linked to a protein
    
    max_linked_unique_per_psm = 0
    max_linked_shared_per_psm = 0
    max_linked_balanced_unique_per_psm = 0
    max_linked_balanced_shared_per_psm = 0
    
    
    for oPsm in lPsm:
        num_unique_per_psm = 0
        num_shared_per_psm = 0
        num_balanced_shared_per_psm = 0
        max_per_psm = 0
        max_balanced_per_psm = 0
        max_unique_per_psm = 0
        max_shared_per_psm = 0
        max_balanced_shared_per_psm = 0
        # get the protein list first
        if oPsm.OriginalPeptide not in original_pep_pro_dict:
            print 'error peptide not in the dictionary'
        else:
            pro_list = original_pep_pro_dict[oPsm.OriginalPeptide]
            for pro_str in pro_list:
                # debug
                num_unique_per_pro = 0
                num_shared_per_pro = 0
                num_balanced_shared_per_pro = 0
                num_per_pro = 0
                num_balanced_per_pro = 0
                # debug
                if pro_str not in pro_original_pep_dict:
                    print 'error in number of sibling peptide'
                else:
                    for pep_str in pro_original_pep_dict[pro_str]:
                        if pep_str not in original_pep_pro_dict:
                            print 'error in pep dictionary'
                        else:
                            if len(original_pep_pro_dict[pep_str]) == 1: # unique
                                num_unique_per_pro += 1
                            else: # shared
                                num_shared_per_pro += 1
                                num_balanced_shared_per_pro += 1.0/float(len(original_pep_pro_dict[pep_str]))
                    num_per_pro += num_unique_per_pro + num_shared_per_pro
                    num_balanced_per_pro += num_unique_per_pro + num_balanced_shared_per_pro
                num_unique_per_psm += num_unique_per_pro
                num_shared_per_psm += num_shared_per_pro
                num_balanced_shared_per_psm += num_balanced_shared_per_pro
                if num_per_pro > max_per_psm:
                    max_per_psm = num_per_pro
                    max_linked_unique_per_psm = num_unique_per_pro
                    max_linked_shared_per_psm = num_shared_per_pro
                if num_balanced_per_pro > max_balanced_per_psm:
                    max_balanced_per_psm = num_balanced_per_pro
                    max_linked_balanced_unique_per_psm = num_unique_per_pro
                    max_linked_balanced_shared_per_psm = num_balanced_shared_per_pro
                if num_unique_per_pro > max_unique_per_psm:
                    max_unique_per_psm = num_unique_per_pro
                if num_shared_per_pro > max_shared_per_psm:
                    max_shared_per_psm = num_shared_per_pro
                if num_balanced_shared_per_pro > max_balanced_shared_per_psm:
                    max_balanced_shared_per_psm = num_balanced_shared_per_pro
            
            if len(original_pep_pro_dict[oPsm.OriginalPeptide]) == 1:
                num_unique_per_psm -= 1
            else:
                num_shared_per_psm -= 1
        
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_shared_per_psm
        oPsm.SPSC = num_shared_per_psm + num_unique_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_shared_per_psm
        oPsm.SPSC = max_shared_per_psm + max_unique_per_psm
        oPsm.SPSC = max_per_psm
        oPsm.UPSC = max_linked_unique_per_psm
        oPsm.SPSC = max_linked_shared_per_psm
        oPsm.SPSC = max_linked_shared_per_psm + max_linked_unique_per_psm
        # balanced
        # oPsm.UPSC = num_unique_per_psm
        # oPsm.SPSC = num_balanced_shared_per_psm
        # oPsm.SPSC = num_balanced_shared_per_psm + num_unique_per_psm
        oPsm.UPSC = max_linked_balanced_unique_per_psm
        oPsm.SPSC = max_linked_balanced_shared_per_psm
        oPsm.SPSC = max_linked_balanced_unique_per_psm + max_linked_balanced_shared_per_psm
        oPsm.SPSC = max_balanced_per_psm
        oPsm.SPSC = max_balanced_shared_per_psm + max_unique_per_psm
    
    

# generate OPSC (# sibling modification), IPSC (# sibling ions, charge), UPSC (# sibling peptides), NMC (# missed cleavage sites)
def generate_Prophet_features(lPsm, config_dict):
    # peptide with PTM dictionary is for IPSC
    peptide_with_modification_dict = {}
    # peptide without PTM dictionary is for OPSC
    peptide_dict = {}
    peptide_protein_dict = {}
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[cleave_after_residues_str],
                                                 config_dict[cleave_before_residues_str])
        if oPsm.IdentifiedPeptide in peptide_with_modification_dict:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] += 1
        else:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] = 1
        if oPsm.OriginalPeptide in peptide_protein_dict:
            pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
            for protein in oPsm.protein_list:
                if not protein in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            pro_list.extend(oPsm.protein_list)
            peptide_protein_dict[oPsm.OriginalPeptide] = pro_list
            

    pattern = re.compile('[^\w\[\]]')
    for key, _value in peptide_with_modification_dict.iteritems():
        peptide_str = pattern.sub('', key)
        if peptide_str in peptide_dict:
            peptide_dict[peptide_str] += 1
        else:
            peptide_dict[peptide_str] = 1
    
    # # sibling peptides
    pro_unique_dict = {}
    pro_shared_dict = {}
    # debug
    pro_balanced_shared_dict = {}
    # debug
    changed_flag = False
    for oPsm in lPsm:
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        changed_flag = False
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                oPsm.protein_list.append(protein)
                changed_flag = True
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
        if len(pro_list) > 1:
            num_pro_f = float(len(pro_list))
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
                # debug
                if protein in pro_balanced_shared_dict:
                    pro_balanced_shared_dict[protein] += 1.0/num_pro_f
                else:
                    pro_balanced_shared_dict[protein] = 1.0/num_pro_f
                # debug
        else:
            if pro_list[0] in pro_unique_dict:
                pro_unique_dict[pro_list[0]] += 1
            else:
                pro_unique_dict[pro_list[0]] = 1
    
    # collect features
    for oPsm in lPsm:
        oPsm.IPSC = peptide_with_modification_dict[oPsm.IdentifiedPeptide] - 1
        oPsm.OPSC = peptide_dict[oPsm.OriginalPeptide] - 1
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        for protein in pro_list:
            if protein in pro_unique_dict:
                oPsm.UPSC += pro_unique_dict[protein]
                if len(pro_list) == 1:
                    oPsm.UPSC -= 1
                
            if protein in pro_shared_dict:
                oPsm.SPSC += pro_shared_dict[protein]
                if len(pro_list) > 1:
                    oPsm.SPSC -= 1
                    
def generate_Prophet_features_group(psm_list, config_dict):
    pep_spectra_U_dict = {}
    pep_spectra_M_dict = {}
    pep_spectra_D_dict = {}
    # iLocalRank = 0, 1, 2, 3
    pep_spectra_UMD_list = [pep_spectra_U_dict, pep_spectra_M_dict, pep_spectra_D_dict, pep_spectra_D_dict]
    for oPsm in psm_list:
        if oPsm.OriginalPeptide not in pep_spectra_UMD_list[oPsm.iLocalRank]:
            pep_spectra_UMD_list[oPsm.iLocalRank][oPsm.OriginalPeptide] = 1
        else:
            pep_spectra_UMD_list[oPsm.iLocalRank][oPsm.OriginalPeptide] += 1
    pro_spectra_U_dict = {}
    pro_spectra_M_dict = {}
    pro_spectra_D_dict = {}
    pro_spectra_UMD_list = [pro_spectra_U_dict, pro_spectra_M_dict, pro_spectra_D_dict, pro_spectra_D_dict]
    for oPsm in psm_list:
        pro_list = oPsm.protein_list
        for pro in pro_list:
            if pro not in pro_spectra_UMD_list[oPsm.iLocalRank]:
                pro_spectra_UMD_list[oPsm.iLocalRank][pro] = 1
            else:
                pro_spectra_UMD_list[oPsm.iLocalRank][pro] += 1
    
    for oPsm in psm_list:
        for i in range(3):
            if oPsm.OriginalPeptide in pep_spectra_UMD_list[i]:
                # if oPsm.OPSC_UMD[i] < pep_spectra_UMD_list[i][oPsm.OriginalPeptide]:
                oPsm.OPSC_UMD[i] = pep_spectra_UMD_list[i][oPsm.OriginalPeptide]
            
        pro_list = oPsm.protein_list

        for pro in pro_list:
            for i in range(3):
                if pro in pro_spectra_UMD_list[i]:
                    if oPsm.SPSC_UMD[i] < pro_spectra_UMD_list[i][pro]:
                        oPsm.SPSC_UMD[i] = pro_spectra_UMD_list[i][pro]



def generate_Prophet_features_test(lPsm, config_dict):
    # peptide with PTM dictionary is for IPSC
    peptide_with_modification_dict = {}
    # peptide without PTM dictionary is for OPSC
    peptide_dict = {}
    peptide_protein_dict = {}
    # psm_set = Set()
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[cleave_after_residues_str],
                                                 config_dict[cleave_before_residues_str])
        '''
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_set:
            continue
        else:
            psm_set.add(unique_id_str)
        '''            
        if oPsm.IdentifiedPeptide in peptide_with_modification_dict:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] += 1
        else:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] = 1
        if oPsm.OriginalPeptide in peptide_protein_dict:
            pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
            for protein in oPsm.protein_list:
                if not protein in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            for pro in oPsm.protein_list:
                if pro not in pro_list:
                    pro_list.extend(oPsm.protein_list)
            peptide_protein_dict[oPsm.OriginalPeptide] = pro_list
            

    pattern = re.compile('[^\w\[\]]')
    for key, _value in peptide_with_modification_dict.iteritems():
        peptide_str = pattern.sub('', key)
        if peptide_str in peptide_dict:
            # peptide_dict[peptide_str] += 1
            peptide_dict[peptide_str] += _value
        else:
            # peptide_dict[peptide_str] = 1
            peptide_dict[peptide_str] = _value
    
    # # sibling peptides
    pro_unique_dict = {} # number of unique peptide of a protain
    pro_shared_dict = {} # number of shared peptide of a protain
    # debug
    pro_balanced_shared_dict = {}
    # debug
    num_changed = 0
    changed_flag = False
    # psm_set.clear()
    pro_pep_dict = {}
    pro_unique_pep_dict = {}
    pro_shared_pep_dict = {}
    for oPsm in lPsm:
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        changed_flag = False
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                oPsm.protein_list.append(protein)
                changed_flag = True
        if len(oPsm.protein_list) != len(pro_list):
            print 'check 3'
        
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
            num_changed += 1
        '''
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_set:
            continue
        else:
            psm_set.add(unique_id_str)
        '''
        if len(pro_list) > 1:
            num_pro_float = float(len(pro_list))
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
                # debug
                if protein in pro_balanced_shared_dict:
                    pro_balanced_shared_dict[protein] += 1.0/num_pro_float
                else:
                    pro_balanced_shared_dict[protein] = 1.0/num_pro_float
                # debug
        else:
            if pro_list[0] in pro_unique_dict:
                pro_unique_dict[pro_list[0]] += 1
            else:
                pro_unique_dict[pro_list[0]] = 1
        
        for pro in pro_list:
            if pro in pro_pep_dict:
                l = pro_pep_dict[pro]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_pep_dict[pro] = l
        if len(pro_list) == 1:
            if pro_list[0] in pro_unique_pep_dict:
                l = pro_unique_pep_dict[pro_list[0]]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_unique_pep_dict[pro_list[0]] = l
        else:
            for pro in pro_list:
                if pro in pro_shared_pep_dict:
                    l = pro_shared_pep_dict[pro]
                    if oPsm.OriginalPeptide not in l:
                        l.append(oPsm.OriginalPeptide)
                else:
                    l = []
                    l.append(oPsm.OriginalPeptide)
                    pro_shared_pep_dict[pro] = l
    print "num changed %d" % num_changed
    
    # collect features
    num_unique_per_pro = 0
    num_shared_per_pro = 0
    num_balanced_shared_per_pro = 0
    
    max_unique_per_psm = 0
    max_shared_per_psm = 0
    max_balanced_shared_per_psm = 0
    
    num_unique_per_psm = 0
    num_shared_per_psm = 0
    num_balanced_shared_per_psm = 0
    
    max_per_psm = 0 # linked to protein
    max_balanced_per_psm = 0 # linked to a protein
    
    max_linked_unique_per_psm = 0
    max_linked_shared_per_psm = 0
    max_linked_balanced_unique_per_psm = 0
    max_linked_balanced_shared_per_psm = 0
    
    for oPsm in lPsm:
        oPsm.IPSC = peptide_with_modification_dict[oPsm.IdentifiedPeptide]
        oPsm.OPSC = peptide_dict[oPsm.OriginalPeptide]
        
        max_unique_per_psm = 0
        max_shared_per_psm = 0
        max_balanced_shared_per_psm = 0
    
        num_unique_per_psm = 0
        num_shared_per_psm = 0
        num_balanced_shared_per_psm = 0
        
        max_per_psm = 0 # linked to protein
        max_balanced_per_psm = 0 # linked to a protein
    
        max_linked_unique_per_psm = 0
        max_linked_shared_per_psm = 0
        max_linked_balanced_unique_per_psm = 0
        max_linked_balanced_shared_per_psm = 0
        
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        
        for protein in pro_list:
            if len(pro_pep_dict[protein]) > oPsm.PPC:
                oPsm.PPC = len(pro_pep_dict[protein])
            if len(pro_list) == 1 and len(pro_unique_pep_dict[protein]) > oPsm.UPPC:
                oPsm.UPPC = len(pro_unique_pep_dict[protein])
            if len(pro_list) > 1 and len(pro_shared_pep_dict[protein]) > oPsm.SPPC:
                oPsm.SPPC = len(pro_shared_pep_dict[protein])
            num_unique_per_pro = 0
            num_shared_per_pro = 0
            num_balanced_shared_per_pro = 0
            
            if protein in pro_unique_dict:
                num_unique_per_pro = pro_unique_dict[protein]
                if len(pro_list) == 1:
                    # pass
                    num_unique_per_pro -= 1
                
            if protein in pro_shared_dict:
                num_shared_per_pro = pro_shared_dict[protein]
                if len(pro_list) > 1:
                    # pass
                    num_shared_per_pro -= 1
            
            if protein in pro_balanced_shared_dict:
                num_balanced_shared_per_pro = pro_balanced_shared_dict[protein]
                if len(pro_list) > 1:
                    num_balanced_shared_per_pro -= 1.0/ float(len(pro_list))
            
            if num_unique_per_pro > max_unique_per_psm:
                max_unique_per_psm = num_unique_per_pro
            if num_shared_per_pro > max_shared_per_psm:
                max_shared_per_psm = num_shared_per_pro
            if num_unique_per_pro + num_shared_per_pro > max_per_psm:
                max_per_psm = num_unique_per_pro + num_shared_per_pro
                max_linked_unique_per_psm = num_unique_per_pro
                max_linked_shared_per_psm = num_shared_per_pro
            num_unique_per_psm += num_unique_per_pro
            num_shared_per_psm += num_shared_per_pro
            
            if num_balanced_shared_per_pro > max_balanced_shared_per_psm:
                max_balanced_shared_per_psm = num_balanced_shared_per_pro
            if num_unique_per_pro + num_balanced_shared_per_pro > max_balanced_per_psm:
                max_balanced_per_psm = num_unique_per_pro + num_balanced_shared_per_pro
                max_linked_balanced_unique_per_psm = num_unique_per_pro
                max_linked_balanced_shared_per_psm = num_balanced_shared_per_pro
            num_balanced_shared_per_psm += num_balanced_shared_per_pro
        
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_shared_per_psm
        oPsm.SPSC = num_unique_per_psm + num_shared_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_shared_per_psm
        # oPsm.SPSC = max_shared_per_psm + max_unique_per_psm
        # oPsm.PPC = oPsm.UPPC + oPsm.SPPC
        # oPsm.UPSC = max_linked_unique_per_psm
        # oPsm.SPSC = max_linked_shared_per_psm
        # oPsm.SPSC = max_linked_balanced_shared_per_psm
        if len(oPsm.protein_list) == 1:
            oPsm.UPSC = 1
        else:
            oPsm.UPSC = 0
        oPsm.SPSC = max_linked_unique_per_psm + max_linked_shared_per_psm
        # oPsm.SPSC = float(max_linked_unique_per_psm)/float(len(oPsm.protein_list)) + float(max_linked_shared_per_psm)/float(len(oPsm.protein_list))
        # oPsm.SPSC = max_per_psm
        
        # oPsm.UPSC = max_unique_per_psm
        # oPsm.SPSC = max_shared_per_psm
        
        '''
        # # balanced
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_balanced_shared_per_psm
        oPsm.SPSC = num_unique_per_psm + num_balanced_shared_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_balanced_shared_per_psm
        oPsm.UPSC = max_linked_balanced_unique_per_psm
        oPsm.SPSC = max_linked_balanced_shared_per_psm
        # oPsm.SPSC = max_unique_per_psm + max_balanced_shared_per_psm
        # oPsm.SPSC = max_linked_balanced_unique_per_psm + max_linked_balanced_shared_per_psm
        '''
# # Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)

# # Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as _f: pass
    except IOError as _e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")

# defaul value
decoy_prefix = 'Rev_'
min_peptide_per_protein = 2
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'

pep_iden_str = '[Peptide_Identification]'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
FDR_Filtering_str = 'FDR_Filtering'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'
cleave_after_residues_str = 'Cleave_After_Residues'
cleave_before_residues_str = 'Cleave_Before_Residues'

## Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    # Save all config values to dictionary
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)

    # only save protein_identification config info to config_dict
    config_dict[decoy_prefix_str] = decoy_prefix
    config_dict[min_peptide_per_protein_str] = min_peptide_per_protein
    config_dict[min_unique_peptide_per_protein_str] = min_unique_peptide_per_protein
    config_dict[remove_decoy_identification_str] = remove_decoy_identification
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + fasta_database_str):
            config_dict[fasta_database_str] = value
        elif key == (pro_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        elif key == (pro_iden_str + min_peptide_per_protein_str):
            config_dict[min_peptide_per_protein_str] = value
        elif key == (pro_iden_str + min_unique_peptide_per_protein_str):
            config_dict[min_unique_peptide_per_protein_str] = value
        elif key == (pro_iden_str + remove_decoy_identification_str):
            config_dict[remove_decoy_identification_str] = value
        elif key == (pep_iden_str + cleave_after_residues_str):
            config_dict[cleave_after_residues_str] = value
        elif key == (pep_iden_str + cleave_before_residues_str):
            config_dict[cleave_before_residues_str] = value
        elif key == (pro_iden_str + FDR_Filtering_str):
            config_dict[FDR_Filtering_str] = value
        else:
            continue

    # return config dictionary
    return config_dict


class Peptide:
    
    def __init__(self):
        self.IdentifiedPeptide = ''
        self.ParentCharge = ''
        self.OriginalPeptide = ''
        self.ProteinNames = []
        self.ProteinCount = 0
        self.SpectralCount = 0
        self.BestScore = 0.0
        self.PSMs = []
        self.ScanType = []
        self.SearchName = []
        
    def add(self, oPsm):
        self.SpectralCount += 1
        if self.BestScore < oPsm.fPredictProbability:
            self.BestScore = oPsm.fPredictProbability
        self.PSMs.append(oPsm.FileName+'['+str(oPsm.ScanNumber) +']')
        self.ScanType.append(oPsm.ScanType)
        self.SearchName.append(oPsm.SearchName)
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        
    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide
        self.ProteinNames = oPsm.ProteinNames
        self.ProteinCount = len(oPsm.protein_list)
        self.SpectralCount = 1
        self.BestScore = oPsm.fPredictProbability
        self.PSMs.append(oPsm.FileName+'['+str(oPsm.ScanNumber) +']')
        self.ScanType.append(oPsm.ScanType)
        self.SearchName.append(oPsm.SearchName)
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'
        
    def __repr__(self):
        l = [self.IdentifiedPeptide,
             str(self.ParentCharge),
             self.OriginalPeptide,
             self.ProteinNames,
             str(self.ProteinCount),
             self.TargetMatch,
             str(self.SpectralCount),
             str(self.BestScore),
             ('{'+','.join(self.PSMs)+'}'),
             ('{'+','.join(self.ScanType)+'}'),
             ('{'+','.join(self.SearchName)+'}')]
        
        return '\t'.join(l) 


def generate_psm_pep_txt(input_file, out_folder, psm_filtered_list):
    # get the FDR, # target, # decoy
    psm_target_int = 0
    psm_decoy_int = 0
    pep_target_int = 0
    pep_decoy_int = 0
    pep_set = Set()
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRev:
            continue
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in pep_set:
            if oPsm.RealLabel == LabelFwd:
                pep_target_int += 1
            else:
                pep_decoy_int += 1
            pep_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            psm_target_int += 1
        else:
            psm_decoy_int += 1
    
    # write out into files
    base_out_filename = input_file.split('/')[-1]
    base_out = out_folder + base_out_filename
    with open(base_out + ".psm.txt", 'w') as fw:
        
        fw.write('#    ########################################\n')
        fw.write('#    ####### PSM Filtering by Sipros ########\n')
        fw.write('#    ########################################\n')
        fw.write('#    \n')
        fw.write('#    #######################\n')
        fw.write('#    # Statistical Results #\n')
        fw.write('#    #######################\n')
        fw.write('#    \n')
        fw.write('#    [Statistical_Results]\n')
        fw.write('#    \n')
        fw.write('#    # Numbers of psm after filtering\n')
        fw.write('#    Decoy_PSMs_After_Filtering = %d\n' % psm_decoy_int)
        fw.write('#    Target_PSMs_After_Filtering = %d\n' % psm_target_int)
        fw.write('#    # PSM FDR = Decoy_PSMs_After_Filtering / Target_PSMs_After_Filtering\n')
        fw.write('#    PSM_FDR = %.2f%%\n' % (100.0*psm_decoy_int/psm_target_int))
        fw.write('#    \n')
        # for psm out
        psm_out_list = ['Filename',  # 0
                    'ScanNumber',  # 1
                    'ParentCharge',  # 2
                    'MeasuredParentMass',  # 3
                    'CalculatedParentMass',  # 4
                    'MassErrorDa',  # 5 CalculatedParentMass - MeasuredParentMass
                    'MassErrorPPM',  # 6 MassErrorDa / CalculatedParentMass
                    'ScanType',  # 7
                    'SearchName',  # 8
                    'ScoringFunction',  # 9
                    'Score',  # 10
                    'DeltaZ',  # 11 the difference score between the rank 1 and 2
                    'DeltaP',  # 12
                    'IdentifiedPeptide',  # 13
                    'OriginalPeptide',  # 14
                    'ProteinNames',  # 15
                    'ProteinCount',  # 16
                    'TargetMatch']  # 17
        fw.write('\t'.join(psm_out_list) + '\n')
        for oPsm in psm_filtered_list:
            if oPsm.RealLabel == LabelRev:
                continue 
            oPsm.clean_protein_name()
            fw.write(oPsm.FileName)
            fw.write('\t')
            fw.write(str(oPsm.ScanNumber))
            fw.write('\t')
            fw.write(str(oPsm.ParentCharge))
            fw.write('\t')
            fw.write('%.3f' % oPsm.MeasuredParentMass)
            fw.write('\t')
            fw.write('%.3f' % oPsm.CalculatedParentMass)
            fw.write('\t')
            fw.write('%.3f' % (oPsm.fMassDiff))
            fw.write('\t')
            fw.write('%.3f' % (1000000*(oPsm.fMassDiff)/oPsm.CalculatedParentMass))
            fw.write('\t')
            fw.write(oPsm.ScanType)
            fw.write('\t')
            fw.write(oPsm.SearchName)
            fw.write('\t')
            fw.write('Sipros10')
            fw.write('\t')
            fw.write('NA')
            fw.write('\t')
            fw.write('NA')
            fw.write('\t')
            fw.write(oPsm.DeltaP)
            fw.write('\t')
            fw.write(oPsm.IdentifiedPeptide)
            fw.write('\t')
            fw.write(oPsm.OriginalPeptide)
            fw.write('\t')
            fw.write(oPsm.ProteinNames)
            fw.write('\t')
            fw.write(str(len(oPsm.protein_list)))
            fw.write('\t')
            if oPsm.RealLabel == LabelFwd:
                fw.write('T')
            else:
                fw.write('F')
            fw.write('\n')
            
    # pep_sub_dict for preparing pep_out
    pep_sub_dict = {}    # initialize dict of list
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRev:
            continue
        pep_ID = oPsm.IdentifiedPeptide + '_+_' + str(oPsm.ParentCharge)
        if pep_ID in pep_sub_dict:
            pep_sub_dict[pep_ID].add(oPsm)
        else:
            oPeptide = Peptide()
            oPeptide.set(oPsm)
            pep_sub_dict[pep_ID] = oPeptide
    
    with open(base_out + ".pep.txt", 'w') as fw:
        # statistic results
        fw.write('#    ########################################\n')
        fw.write('#    ####### PSM Filtering by Sipros ########\n')
        fw.write('#    ########################################\n')
        fw.write('#    \n')
        fw.write('#    #######################\n')
        fw.write('#    # Statistical Results #\n')
        fw.write('#    #######################\n')
        fw.write('#    \n')
        fw.write('#    [Statistical_Results]\n')
        fw.write('#    \n')
        fw.write('#    # Numbers of peptide after filtering\n')
        fw.write('#    Decoy_peptides_After_Filtering = %d\n' % pep_decoy_int)
        fw.write('#    Target_peptides_After_Filtering = %d\n' % pep_target_int)
        fw.write('#    # peptide FDR = Decoy_peptides_After_Filtering / Target_peptides_After_Filtering\n')
        fw.write('#    peptide_FDR = %.2f%%\n' % (100.0*pep_decoy_int/pep_target_int))
        fw.write('#    \n')
        
        # for pep out
        pep_out_list = ['IdentifiedPeptide',    #0
                    'ParentCharge',         #1
                    'OriginalPeptide',      #2
                    'ProteinNames',         #3
                    'ProteinCount',         #4
                    'TargetMatch',          #5
                    'SpectralCount',        #6 number of PSMs matched to this peptide
                    'BestScore',            #7 the highest score of those PSMs
                    'PSMs',                 #8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                    'ScanType',             #9
                    'SearchName']           #10
        fw.write('\t'.join(pep_out_list) + '\n')
        for _pep_id, oPeptide in pep_sub_dict.iteritems():
            fw.write(repr(oPeptide))
            fw.write('\n')


## check pep and psm files pair set, and save run#
def get_run_num(pep_file_list, psm_file_list):

    # dictionary of Run# for each pep_file
    run_num_dict = {}
    psm_run_num_dict = {}

    # read multiple pep files
    for pep_file_idx, pep_file in enumerate(pep_file_list):
    
        # check file exist and open
        check_file_exist(pep_file)

        # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
        base_pep_file = pep_file.replace(pep_file_ext, "")

        # make a psm filename using pep filename
        psm_file = base_pep_file + psm_file_ext

        # check psm file exist
        check_file_exist(psm_file)

        # for Run#
        run_tag = 'Run' + str(pep_file_idx + 1)

        # run_num_dict
        run_num_dict[pep_file] = run_tag
        psm_run_num_dict[psm_file] = run_tag

    return (run_num_dict, psm_run_num_dict)

PepOutFields = sipros_post_module.PepOutFields
# # need pep.txt file and pro.txt file, and big PSM table file
# # Spectral Count for original/identified peptide
# # Unique peptide counts and total peptide counts for protein
def feature_update(run_num_dict, pro_file, psm_tab, psm_tab_new):
    
    # save the pep file and pro file data to the defaultdict
    id_pep_data_dict = {}
    or_pep_data_dict = {}
    pro_data_dict = {}
    
    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, _run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
        f = open(pep_file, 'rb')
        # read line with csv
        pep_reader = csv.reader(CommentedFile(f),
                                   delimiter='\t')
        # skip header
        _headline = pep_reader.next()
        
        # get data
        for pep_line in pep_reader:
            # adapt class PepOutFields
            pep_obj = PepOutFields._make(pep_line)
            identified_peptide = pep_obj.IdentifiedPeptide
            identified_peptide = identified_peptide.strip()
            ParentCharge = pep_obj.ParentCharge
            identified_peptide = identified_peptide +'_' + ParentCharge
            original_peptide = pep_obj.OriginalPeptide
            original_peptide = original_peptide.strip()
            
            iSpectralCount = int(pep_obj.SpectralCount.strip())
            
            if identified_peptide in id_pep_data_dict:
                id_pep_data_dict[identified_peptide] += iSpectralCount
            else:
                id_pep_data_dict[identified_peptide] = iSpectralCount
            
            if original_peptide in or_pep_data_dict:
                or_pep_data_dict[original_peptide] += iSpectralCount
            else:
                or_pep_data_dict[original_peptide] = iSpectralCount
        f.close()
    
    f = open(pro_file, 'rb')
    pro_reader = csv.reader(CommentedFile(f),
                                   delimiter='\t')
    # skip header
    headline = pro_reader.next()
    iNumRuns = (len(headline) - 2)/6
    iUniquePeptideCounts = 0
    iTotalPeptideCounts = 0
    for asWords in pro_reader:
        identified_protain = asWords[0]
        iUniquePeptideCounts = 0
        iTotalPeptideCounts = 0
        for i in range(iNumRuns):
            iUniquePeptideCounts += int(asWords[i*6+1])
            iTotalPeptideCounts += int(asWords[i*6+2])
        pro_data_dict[identified_protain] = [iUniquePeptideCounts, iTotalPeptideCounts]
    
    f.close()
    
    fr = open(psm_tab, 'rb')
    psm_reader = csv.reader(CommentedFile(fr),
                                   delimiter='\t')
    
    # skip header
    headline = psm_reader.next()
    
    with open(psm_tab_new, 'w') as f:
        # header
        f.write('\t'.join(headline))
        f.write('\tpep_psm\t')
        f.write('pro_pep')
        f.write('\n')
        
        for asWords in psm_reader:
            original_peptide = asWords[7]
            if original_peptide in or_pep_data_dict:
                iNumPepPsm = or_pep_data_dict[original_peptide]
                iNumPepPsm -= 1
                if iNumPepPsm < 0:
                    print 'Error'
                    exit(1)
            else:
                iNumPepPsm = 0
            sProteins = asWords[12]
            asProteins = (sProteins[1:-1]).split()
            iNumTotalPep = 0
            iNumUniqPep = 0
            for sProtein in asProteins:
                if sProtein in pro_data_dict:
                    l = pro_data_dict[sProtein]
                    iNumTotalPep += l[1]
                    iNumUniqPep += l[0]
            f.write('\t'.join(asWords))
            f.write('\t')
            f.write(str(iNumPepPsm))
            f.write('\t')
            if iNumUniqPep > 1:
                f.write('2')
            elif iNumTotalPep > 1:
                f.write('1')
            else:
                f.write('0')
            f.write('\n')
    fr.close()
    
def clean_folder(output_folder):
    for the_file in os.listdir(output_folder):
        file_path = os.path.join(output_folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            # elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

## Glboal variables
pep_file_ext = '.pep.txt'
psm_file_ext = '.psm.txt'
read_fasta_file = sipros_peptides_assembling.read_fasta_file
get_file_list_with_ext = sipros_post_module.get_file_list_with_ext
get_base_out = sipros_post_module.get_base_out
read_run_files = sipros_peptides_assembling.read_run_files
greedy_alg = sipros_peptides_assembling.greedy_alg
report_output = sipros_peptides_assembling.report_output

def assembly(output_folder, config_dict, psm_tab_file):
    # Read fasta file and retrieve protein ID and description
    fasta_ID_dict = read_fasta_file(output_folder, config_dict)
    
    # Get .pep.txt output file(s) in working directory
    pep_file_list = get_file_list_with_ext(output_folder, pep_file_ext)
    # Get .psm.txt output file(s) in working directory
    psm_file_list = get_file_list_with_ext(output_folder, psm_file_ext)
    # check pep and psm files pair set, and save run#
    (run_num_dict, psm_run_num_dict) = get_run_num(pep_file_list, psm_file_list)
    # Get base_out for output
    base_out_default = 'Sipros_searches'
    base_out = get_base_out(pep_file_list, base_out_default, output_folder)
    sys.stderr.write('Done!\n')
    
    # Read and load pep and psm files
    sys.stderr.write('[Step 2] Load %s file(s): Running -> ' % (pep_file_ext))
    (pep_data_dict, psm_data_dict, pro_pep_dict, pep_pro_dict) = read_run_files(run_num_dict)
    sys.stderr.write('Done!\n')

    # Merge indistinguishable proteins that have an identical set of peptides
    sys.stderr.write('[Step 3] Merge indistinguishable proteins: Running -> Done!\n')

    # extract proteins that have >2 peptides and at least one of those is unique
    # then iteratively extract a protein at a time that covers the most peptides
    sys.stderr.write('[Step 4] Greedy algorithm for identifying a list of proteins: Running -> ')
    (pro_greedy_list) = greedy_alg(config_dict, pro_pep_dict, pep_pro_dict)
    sys.stderr.write('Done!\n')
    
    # Report output
    sys.stderr.write('[Step 5] Report output: Running -> ')

    # Report output files
    report_output(config_dict,
                  run_num_dict,
                  psm_run_num_dict,
                  pep_data_dict,
                  psm_data_dict,
                  pro_pep_dict,
                  pep_pro_dict,
                  pro_greedy_list,
                  fasta_ID_dict,
                  base_out)
    sys.stderr.write('Done!\n')
    return None
    # Feature Update
    sys.stderr.write('[Step 6] Feature update: Running -> ')
    base_out_filename = psm_tab_file.split('/')[-1]
    psm_tab_new = output_folder + base_out_filename
    feature_update(run_num_dict, base_out+'.pro.txt', psm_tab_file, psm_tab_new)
    return psm_tab_new

def show_measured_predicted_rt(psm_list, filename_prefix):
    fw_fwr = open(filename_prefix+"_fwr.txt", 'w')
    fw_fwr.write("measuread\tpredicted\n")
    fw_shu = open(filename_prefix+"_shu.txt", 'w')
    fw_shu.write("measuread\tpredicted\n")
    
    for oPsm in psm_list:
        if oPsm.RealLabel == LabelFwd:
            fw_fwr.write(str(oPsm.fRtMeasured))
            fw_fwr.write('\t')
            fw_fwr.write(str(oPsm.fRtPredict))
            fw_fwr.write('\n')
        else:
            fw_shu.write(str(oPsm.fRtMeasured))
            fw_shu.write('\t')
            fw_shu.write(str(oPsm.fRtPredict))
            fw_shu.write('\n')
    
    fw_fwr.close()
    fw_shu.close()
    
def stacking(block_num_int, psm_list):
    # # construct training data
    psm_filtered_list = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []     
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\n")
    label_np = np.array(label_list)
    
    # split the data into B blocks
    train_data_block = [[], []]
    train_label_block = [[], []]
    block_idx_list = np.random.randint(0, block_num_int, size=len(data_np))
    data_idx_int = 0
    for c in block_idx_list:
        if c == 0:
            train_data_block[0].append(data_np[data_idx_int])
            train_label_block[0].append(label_np[data_idx_int])
        else:
            train_data_block[1].append(data_np[data_idx_int])
            train_label_block[1].append(label_np[data_idx_int])
        data_idx_int += 1
        
    # features and models
    feature_matrix = [[1, 2, 3, 15, 16, 17],    # score section
                      [5, 29],                  # mass section
                      [24],                     # digestion section
                      [25, 26, 27, 28],         # pep pro support section
                      [30]]                     # PTM section
    all_feature_list = [1, 2, 3, 5, 15, 16, 17, 24, 25, 26, 27, 28, 29, 30]
    tier_1 = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        # tier_1[i].fit(np.array(train_data_block[1])[:, feature_matrix[i]], np.array(train_label_block[1]))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
        predict_np = tier_1[i].predict(data_np[:, feature_matrix[i]])
        # show_TP_TN_FP_FN(label_np, predict_np)
    
    
    
    
    tier_2 = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    # # train tier_2
    train_tier_2_np = []
    for i in range(len(feature_matrix)):
        predict_np = tier_1[i].predict_proba(data_np[:, feature_matrix[i]])
        train_tier_2_np.append(predict_np[:, 1])
    train_tier_2_np = np.transpose(np.array(train_tier_2_np))
    tier_2.fit(train_tier_2_np, label_np)
    '''
    
    # # re-train tier_1
    tier_1 = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
    '''
    # # testing
    unknown_np = np.array(unknown_list)
    test_tier_2_np = []
    for i in range(len(feature_matrix)):
        predict_np = tier_1[i].predict_proba(unknown_np[:, feature_matrix[i]])
        test_tier_2_np.append(predict_np[:, 1])
    test_tier_2_np = np.transpose(np.array(test_tier_2_np))
    predict_np = tier_2.predict_proba(test_tier_2_np)
    
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
    
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    psm_filtered_list_local = show_Fdr(psm_list, None, fdr=0.01)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
        
    return psm_filtered_list

def remark_concensus(psm_list):
    psm_dict = {}
    for oPsm in psm_list:
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_dict:
            psm_dict[unique_id_str] += 1
        else:
            psm_dict[unique_id_str] = 1
    
    psm_set = Set()
    for oPsm in psm_list:
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        count_int = psm_dict[unique_id_str]
        oPsm.iLocalRank = 3 - count_int
        if count_int == 2:
            psm_set.add(oPsm.FileName + '_' + str(oPsm.ScanNumber))
    
    for oPsm in psm_list:
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber)
        if unique_id_str in psm_set:
            if oPsm.iLocalRank == 2:
                oPsm.iLocalRank = 3 # Mi

def mass_filter(psm_list):
    psm_new_list = []
    for oPsm in psm_list:
        if abs(oPsm.fMassDiff) > mass_tolerance:
            continue
        psm_new_list.append(oPsm)
    return psm_new_list

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse options
    (input_file, config_file, output_folder) = parse_options(argv)
    
    # get the configuration parameters
    config_dict = parse_config(config_file)
    
    # read the big psm table
    psm_list = read_psm_table(input_file)
    
    # get score agreement info
    remark_concensus(psm_list)
    
    # mass filtering
    psm_list = mass_filter(psm_list)
    
    # generate features
    # generate_Prophet_features(psm_list, config_dict)
    generate_Prophet_features_test(psm_list, config_dict)
    # generate_pep_pro_features(psm_list, config_dict)
    generate_Prophet_features_group(psm_list, config_dict)
    
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_list()
        
    # machine learning
    #test_random_forest(psm_dict, psm_list)
    del feature_selection_list[:]
    feature_selection_list.extend([1, 2, 3, 5, 15, 16, 17, 24, 26, 28]) #
    
    psm_filtered_list = logistic_regression_no_category(psm_dict, psm_list, 0.01/2, None, output_folder, input_file, config_dict)

    generate_psm_pep_txt(input_file, output_folder, psm_filtered_list)

    print 'Done.'
    
    return
    
def main2(argv=None):
    return
    if argv is None:
        argv = sys.argv

    # parse options
    (input_file, config_file, output_folder, negative_file) = parse_options(argv)
    
    # get the configuration parameters
    config_dict = parse_config(config_file)
    # config_dict[min_peptide_per_protein_str] = 1
    # config_dict[min_unique_peptide_per_protein_str] = 1
    assembly(output_folder, config_dict, input_file)
    print 'Done.'
    return
    

if __name__ == '__main__':
    main()
    sys.exit(main2())
