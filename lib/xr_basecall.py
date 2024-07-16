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

#Generate directories
print('Xemora [Status] - Initializing Xemora basecalling.')
print('??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????')

working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_fastq_dir = check_make_dir(os.path.join(mod_dir,'fastq'))
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
    cmd=os.path.expanduser(basecaller_path)+' -i '+mod_pod_dir+' -s '+mod_fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora  [STATUS] - Skipping POD5 basecalling for modified bases.')
'''
#Step 3: Merge Bam files 
if os.path.isfile(os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam') == False or regenerate_bam == True: 
    cmd = 'samtools merge '+os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam'+' '+os.path.join(mod_fastq_dir,'pass/*.bam -f')
    print('Xemora  [STATUS] - Merging modified BAM files.')
    os.system(cmd)
'''


# Merge modified bam files
if regenerate_bam == True: 
#if os.path.isfile(os.path.join(mod_bam_dir, os.path.basename(mod_bam_dir)) + '.bam') == False or regenerate_bam == True: 

    #Purge directory
    cmd_rmv = 'rm -rf '+mod_bam_dir+'/*'
    os.system(cmd_rmv)
    
    #Set modified bam directory path 
    mod_bam_path = os.path.join(mod_bam_dir, os.path.basename(mod_bam_dir))
    
    # Merging pass bam files
    cmd_pass = 'samtools merge ' + mod_bam_path + '_pass.bam ' + os.path.join(mod_fastq_dir, 'pass/*.bam -f')
    print('Xemora [STATUS] - Merging modified PASS BAM files.')
    os.system(cmd_pass)
    
    print('Amount of total reads in modified PASS BAM file') 
    cmd = 'samtools view -c ' + mod_bam_path + '_pass.bam'
    os.system(cmd) 
    print('Amount of unaligned reads in modified PASS BAM file') 
    cmd = 'samtools view -c -f 4 ' + mod_bam_path + '_pass.bam'
    os.system(cmd)
    
    
    if merge_fail == True:
        # Merging fail bam files
        print('Xemora [STATUS] - Merging modified FAIL BAM files.')
        cmd_fail = 'samtools merge ' + mod_bam_path + '_fail.bam ' + os.path.join(mod_fastq_dir, 'fail/*.bam -f')
        os.system(cmd_fail)
        
        print('Amount of total reads in modified FAIL BAM file')
        cmd = 'samtools view -c ' + mod_bam_path + '_fail.bam'
        os.system(cmd) 
        print('Amount of unaligned reads in modified FAIL BAM file')
        cmd = 'samtools view -c -f 4 ' + mod_bam_path + '_fail.bam'
        os.system(cmd)
        
        
    #Generate merged bam
    cmd_both = 'samtools merge ' + mod_bam_path + '_all.bam ' + mod_bam_path+'*.bam -f'
    os.system(cmd_both)

    print('Amount of total reads in modified FULL BAM file') 
    cmd = 'samtools view -c ' + mod_bam_path + '_all.bam'
    os.system(cmd) 
    print('Amount of unaligned reads in modified FULL BAM file') 
    cmd = 'samtools view -c -f 4 ' + mod_bam_path + '_all.bam'
    os.system(cmd)
else: 
        #Set modified bam directory path 
    mod_bam_path = os.path.join(mod_bam_dir, os.path.basename(mod_bam_dir))
    

#Step 4: Bed file generation 
if os.stat(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))).st_size == 0: 
    print('Xemora  [ERROR] - Empty xfasta file generated. Check that XNA bases were present in sequence of input fasta file.')
    sys.exit()

print('Xemora  [STATUS] - Generating bed file for modified base.')
print(os.path.join(ref_dir,mod_base))

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
            --basecall-anchor \
      --chunk-context '+chunk_context
    os.system(cmd)
    '''
    print('Xemora  [STATUS] - Generating chunks for modified basecalling.')
    #no motif
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5'+' \
      '+mod_bam_path+'_all.bam'+' \
      --output-remora-training-file '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,mod_base)+'.bed'+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --basecall-anchor \
       --motif N 0 \
      --chunk-context '+chunk_context
    os.system(cmd)
try: 
    print('Xemora  [STATUS] - Performing basecalling.')
    cmd = 'remora \
      validate from_remora_dataset \
      '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
      --model '+os.path.expanduser(model_file)+' \
      --full-results-filename '+os.path.join(working_dir,'per-read_modifications.tsv')+' \
      --out-file '+os.path.join(working_dir,'summary_modifications.tsv')
    os.system(cmd)

    print('Xemora  [STATUS] - Basecalling done.')
    print('Xemora  [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'per-read_modifications.tsv'))
    print('Xemora  [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'summary_modifications.tsv'))
    print('Xemora  [STATUS] - Exiting')
except:
    print('Xemora  [ERROR] - Failed to initialize basecalling model. Check logs.')

