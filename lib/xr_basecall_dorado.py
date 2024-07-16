########################################################################
########################################################################
"""
xr_basecall.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
########################################################################

import os
import glob
import sys
from pathlib import Path
from xr_tools  import *
from xr_params import *

############################################################
working_dir = os.path.expanduser(sys.argv[1])
raw_dir = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
model_file = os.path.expanduser(sys.argv[4])
############################################################
###################Parameters##############################
#Method to evaluate modication in dataset (either 'validation' or 'infer')
eval_method = 'validate'

'''
need to edit bed file generation 
'''

#Range of chunk context to use (in bp) for modified base training (default +/- 0) 
mod_chunk_range = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0
#########################################################################
#Generate directories
print('Xemora [Status] - Initializing Xemora basecalling.')

working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
#########################################################################

#Step 0: FASTA to xFASTA conversion
if os.path.isfile(os.path.expanduser(xna_ref_fasta)): 
    cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(xna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora  [ERROR] - Reference fasta xna file not file. Please check file exist or file path.')
    sys.exit()

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

#Step 2: #Basecall pod5 files 
if basecall_pod ==True: 
    #cmd=os.path.expanduser(basecaller_path)+' -i '+mod_pod_dir+' -s '+mod_fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    cmd = '{} basecaller {} --reference {} --no-trim --emit-moves --min-qscore {} {} > {}'.format(dorado_path, dorado_model, os.path.join(ref_dir, 'x'+os.path.basename(xna_ref_fasta)), min_qscore, mod_pod_dir, os.path.join(mod_bam_dir, 'bc.bam'))
    os.system(cmd)
else: 
    print('Xemora  [STATUS] - Skipping POD5 basecalling for modified bases.')

#Step 4: Bed file generation 
if os.stat(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))).st_size == 0: 
    print('Xemora  [ERROR] - Empty xfasta file generated. Check that XNA bases were present in sequence of input fasta file.')
    sys.exit()

print('Xemora  [STATUS] - Generating bed file for modified base.')
print(os.path.join(ref_dir,mod_base))

'''
REMINDER: HARD CODED TO BE +/- 3 BASES RIGHT NOW, MAKE THIS A FUNCTION WITH RANGE AND SHITFTS AS INPUT LATER
'''
cmd = 'python lib/xr_xfasta2bed.py '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))+' '+os.path.join(ref_dir,mod_base+'.bed ' +mod_base+' '+mod_base)
os.system(cmd)

#Step 5: Generate Chunks. 
if regenerate_chunks == True:
    '''
    print('Xemora  [STATUS] - Generating chunks for modified basecalling.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5'+' \
      '+mod_bam_path+'_all.bam'+' \
      --output-remora-training-file '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,mod_base)+'.bed'+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)
    '''
    print('Xemora  [STATUS] - Generating chunks for modified basecalling.')
    #no motif
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5'+' \
      '+os.path.join(mod_bam_dir, 'bc.bam')+' \
      --output-remora-training-file '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,mod_base)+'.bed'+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
       --motif N 0 \
      --chunk-context '+chunk_context
    os.system(cmd)
try: 
    print('Xemora  [STATUS] - Performing basecalling.')
    if eval_method == 'validate':
        cmd = 'remora \
          validate from_remora_dataset \
          '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
          --model '+os.path.expanduser(model_file)+' \
          --full-results-filename '+os.path.join(working_dir,'per-read_modifications.tsv')+' \
          --out-file '+os.path.join(working_dir,'summary_modifications.tsv')
        os.system(cmd)

        print('Xemora [STATUS] - Basecalling done.')
        print('Xemora [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'per-read_modifications.tsv'))
        print('Xemora [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'summary_modifications.tsv'))
        #print('Xemora [STATUS] - Exiting')
        
        
        print('Filtering out per-red_modifications.tsv and generating consensus results')
        #Calculating consensus modification
        cmd = 'python lib/xr_consensus_modifications.py '+working_dir+' '+os.path.join(working_dir,'per-read_modifications.tsv')+' '+os.path.join(ref_dir,mod_base)+'.bed'
        os.system(cmd)

    elif eval_method == 'infer':
        cmd = 'remora \
                infer from_pod5_and_bam \
                --out-bam ' +os.path.join(working_dir, 'results.bam')+' \
                --model '+os.path.expanduser(model_file) +' \
                ' +os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5 \
                ' +os.path.join(mod_bam_dir, 'bc.bam')
        os.system(cmd)
        
        #sam file conversion for manual parsing 
        cmd = 'samtools view -ho '+os.path.join(working_dir,'results.sam') + ' ' + os.path.join(working_dir, 'results.bam')
        os.system(cmd) 
        
        print('Xemora [STATUS] - Basecalling done.')
        print('Xemora [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'results.bam'))
        print('Xemora [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'results.sam'))
        
        
        #Mod kit analysis
        cmd = 'modkit summary --no-sampling --log-filepath '+os.path.join(working_dir,'modkit_logs.csv') + ' ' +os.path.join(working_dir, 'results.bam')
        os.system(cmd)
        
    else:
        print("Xemora [ERROR] - Please set eval_method to 'validate' or 'infer'")
        sys.ext()
except:
    print('Xemora [ERROR] - Failed to initialize basecalling model. Check logs.')

