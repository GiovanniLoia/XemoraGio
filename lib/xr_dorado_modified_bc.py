########################################################################
########################################################################
"""

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
import os
import glob
import sys


from xr_tools  import *
from xr_params import *
########################################################################
#############################Parameters#################################



working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/model_evaluations/dorado_modified_basecalling_tests/trial_with_ST_Model'

raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5'

dorado_model = '/home/marchandlab/dorado-0.7.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0'

#xna_model_path = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/BN_balanced_dorado'
xna_model_path = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/ST_Train_Q9_Dorado'

dorado_path ='/home/marchandlab/dorado-0.7.0-linux-x64/bin/dorado'

min_qscore = 7
########################################################################
#Step 0 - Initialize
print('Xemora [STATUS] - Initializing Xemora...')

##Generate directories
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
mod_pod_dir = check_make_dir(os.path.join(working_dir,'pod5'))
mod_bc_dir = check_make_dir(os.path.join(working_dir,'basecall'))
mod_bam_dir = check_make_dir(os.path.join(working_dir,'bam'))

#######################################################################
'''
#Step 0: FASTA to xFASTA conversion
if os.path.isfile(os.path.expanduser(xna_ref_fasta)): 
    cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(xna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora  [ERROR] - Reference fasta xna file not file. Please check file exist or file path.')
    sys.exit()
'''
#Step 1: Generate or merge pod5 files if needed
xna_file_type = os.listdir(raw_dir)
# Check if the directory is not empty
if xna_file_type:
    # Get the first file in the directory
    xna_first_file = xna_file_type[0]
    # Check if the first file is a .pod5 file
    if xna_first_file.endswith(".pod5"):
        if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            pod5_merge(get_pod5_subdir(raw_dir), os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5')
        else:
            print('Xemora [STATUS] - POD5 files merged. Skipping merge')
    else:
        if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            cod5_to_fast5(get_fast5_subdir(raw_dir), os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5')
        else: 
            print('Xemora [STATUS] - Converted POD5 file for modified base found. Skipping POD5 coversion')
else:
    print('Xemora [ERROR] - Modified Fast5/POD5 directory empty, please check input directory')
    sys.exit()

#dorado command running Xemora model
#cmd = '{} basecaller {} --reference {} --no-trim --emit-moves --min-qscore {} {} > {}'.format(dorado_path, dorado_model, os.path.join(ref_dir, 'x'+os.path.basename(xna_ref_fasta)), min_qscore, mod_pod_dir, os.path.join(mod_bam_dir, 'bc.bam'))
cmd = '{} basecaller {} --modified-bases-models {} --no-trim --emit-moves --min-qscore {} {} > {}'.format(dorado_path, dorado_model, xna_model_path, min_qscore, mod_pod_dir, os.path.join(mod_bam_dir, 'bc.bam'))
os.system(cmd)

#convert to sam for manual parsing 
cmd = 'samtools view -ho '+os.path.join(mod_bam_dir, 'bc.sam')+' '+os.path.join(mod_bam_dir, 'bc.bam')
