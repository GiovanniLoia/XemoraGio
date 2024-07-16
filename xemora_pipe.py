########################################################################
########################################################################
"""
xemora_pipe.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
########################################################################


import os
import numpy as np
import glob
import sys
from pathlib import Path
from lib.xr_tools import *
from lib.xr_params import *


############################################################
#Training paths
'''
BSn/AT 90-mer variations 
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/240326_NCNANGBCNCNTN_alternating_no_validate_edits/S' #Input directory where data processing/outputs will be for training
xna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #Input directory for xna containing fast5 or pod5 files
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa' #Input an xna containing fasta file for modified dataset above, NNNBNNN
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train.fa' #Input an xna containing fasta file for canonical dataset above 
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_changing_randomer.fa' 
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_randomer_ends.fa' #NNN [...] B [...] NNN
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_alternate_N.fa' #...NCNANGBCNCNTN...
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_NNN_constant3_B.fa' #...NNNAGGBCTCNNN...
dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240109_AT_90mer_xr_train_rerun/50f5' #Input directory for canonical analogous containing fast5 or pod5 files 
#dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240109_AT_90mer_xr_train_rerun/fast5_50-72_basecall'
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa' #Input an xna containing fasta file for canonical dataset above, NNNBNNN
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train.fa' #Input an xna containing fasta file for canonical dataset above 
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_changing_randomer.fa' 
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_alternate_N.fa' #...NCNANGBCNCNTN..
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_NNN_constant3_B.fa' #...NNNAGGBCTCNNN...
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_randomer_ends.fa' #NNN [...] B [...] NNN
'''

#PZ 139-mer/GC 131-mer variations
#working_dir =  '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/PZ_variation_tests/240407_NNN_constant3_P_constant3_NNN/Z'
#working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/PZ_variation_tests/240404_null_test/P'
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/PZ_variation_tests/number_of_Ns'
xna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/20240215_1810_MN37138_ARS988_4bbd5246/pod5_0-72'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/reference/PZ_NB25_xr_Train_3N.fasta' #Ground truth
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/reference/PZ_NB25_xr_Train_NNNPNNN.fasta'#NNNPNNN
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/reference/PZ_NB25_xr_Train_randomer_ends.fasta'
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/reference/PZ_NB25_xr_Train_3constant.fasta' #3 bases away from xna 
dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/20240216_1817_MN41475_ASE526_f9fc38c7/100_pod5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/reference/GC_71mer_xr_Train_3N.fasta'
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/reference/GC_71mer_xr_Train_NNNPNNN.fasta'#NNNPNNN
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/reference/GC_71mer_xr_Train_randomer_ends.fasta'
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/reference/GC_71mer_xr_Train_constant3.fasta' #3 bases away from xna
############################################################

############################################################
#Basecall paths

#BSn/AT 90-mer variations 
#bc_working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/model_evaluations/B_single_context_basecaller_model_changes/sup' #Input directory where data processing/outputs will be for basecalling
#bc_working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/model_evaluations/range_tests/B_xenomorph_library_BN_bal_+-5' #Input directory where data processing/outputs will be for basecalling
bc_working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/model_evaluations/range_tests/A_single_context_BN_bal_sup_model_+-30' #Input directory where data processing/outputs will be for basecalling
#bc_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5'

#bc_raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/231003_BSn_libv4_FLG114/20231003_1555_MN37138_AQK018_63a7330b/fast5' #BSn xm lib 
#bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/231003_BSn_libv4_FLG114/ref_libv2_BS_CxDx-.fa' #BSn xm lib (can be used for GC minion run as well

bc_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_AT_90mer_xr_train/20240104_1445_MN41475_ARW614_f45baf6f/fast5' #AT single se
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train.fa' #Input an xna containing fasta file for canonical dataset above 
bc_model_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/BN_balanced_training/BN_model.pt' #Input either provided model files or generated models from xr_train
#bc_model_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/SN_unbalanced_training/model/SN_model.pt' #unbalanced SN model 
#bc_model_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/SN_balanced_training/model/SN_model.pt'
#bc_model_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/BN_Unbalanced/model/BN_model.pt'
############################################################

############################################################
train_model = False
data_vis = False
basecall_reads = True
output_basecall_results = False
############################################################

#conda activate xemora-re
        
#Train dataset with xemora train
if train_model ==True: 
    cmd = 'python xemora.py train -w '+working_dir+' -f '+xna_raw_dir+' '+dna_raw_dir+' -r '+xna_ref_fasta+' '+dna_ref_fasta
    os.system(cmd)

if data_vis ==True:
    #Needs train to be ran before
    cmd = 'python lib/xr_datavis.py '+working_dir+' '+xna_raw_dir+' '+dna_raw_dir+' '+xna_ref_fasta+' '+dna_ref_fasta
    os.system(cmd)
#Basecall fast5 directory 
if basecall_reads==True: 
    print(bc_raw_dir)
    cmd = 'python xemora.py basecall -w '+bc_working_dir+' -f '+bc_raw_dir+' -r '+bc_xna_ref_fasta+' -m '+bc_model_file 
    os.system(cmd)

#output results
if output_basecall_results==True: 
    script_path = './lib/xr_results.py'
    cmd = f'python {script_path} {bc_working_dir}'
    os.system(cmd)
    
