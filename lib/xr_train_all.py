########################################################################
########################################################################
"""
xr_train_all.py 

Title: Unpublished work

By: J. Sumabat, J. A. Marchand

Updated: 5/4/23

Note - fallback directories seem hardcoded.

Change log: 
-Removed calculation of purine/pyrimidine labels. Labels are never used downstream. Not needed. 
    labels also hardcode structural information that is not necessarily relevant 
-Changed variable name "basecall" to "raw_basecall" 
-Changed other variable name "basecall" to "regenerate basecall"
-Rewrote ref_xnapos_gen since it had unused outputs and was hard coding XNA labels
-Rewrote find_xna_position to output XNA label (strand corrected). Removed unecessary lines 
-Fixed a problem in xna_base_pos_all for reverse complement - Needed to rc the raw basecall 
-Filtering now happens before bed file generation and through two steps
>Strand filter (should be obsolete) followed by XNA base label filtering
>XNA base label filtering is required for modified reads
>XNA base label filtering is toggleable for canonical reads
-Can_base and mod_base settings in xr_params have been updated to can_base_train and mod_base_train (training specific)
making it easier to switch between basecalling and trianing 
-Moved execution of pod5 merging/fast5 conversion to merge_convert_reads function, rather than merge_pod_command (deprecated)
-Added mapping score filter to xr_params.py. This is a first alignment filter


-Note - 
-Currently Param files has BN, SN as confounding pairs. 
-Currently set to can_base_train (N) = which means use all standard DNA reads regardless of what position the mask mapped to in first alignment
-Currently set to use all canonical base reads, even if its different context than modified base data. This seems like it is fine and does not affect performance. 


Note on reference file used for training: 
N10 masking works for model training. 
For basecalling, no mask is needed. 

Note on model training: 
Setting sliding range for chunks improves model performance. Model training in this case does not need to be balanced. 
Modified base - If set to 0, this maximizes positional recall
Canonical base - Set to large to use canonical reads for reading regular DNA contexts 

Examples: [Mod range, canonical range]
[0,0] - Standard way we would build models 
[0,4] - Preserves recall of XNA at 0 position, loses accuracy of XNA with increasing range. Generally improves sensitivity as well. 
[4,4] - Improvement in performance across the board for accuracy, at the expense of lower positional resolution. 

Note: [0,4] model or similar might be the best way to train, since it maximizes examples of canonical sets while preserving positional information. 
"""


########################################################################
########################################################################

import os
import glob
import sys
from pathlib import Path
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
from collections import Counter
from xr_tools  import *
from xr_params import *
from xr_all_context_methods import *

############################################################
#extra directories for PZ
'''
#PZ 139 mer
#working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240428_only_G_bed_shift_3_jay'
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240428_PZ_shift_0_no_motif'
xna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/20240215_1810_MN37138_ARS988_4bbd5246/pod5_0-72' #PZ
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/reference/PZ_NB25_xr_Train_3N.fasta' #PZ 
dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/20240216_1817_MN41475_ASE526_f9fc38c7/100_pod5' #PZ
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/reference/GC_71mer_xr_Train_3N.fasta' #PZ
'''
######################Parameters###########################

#toggle basecalling 
regenerate_basecall = True

#If true  only canonical matching mod label base are used (context-matching the same strand). If false all canonical reads used for training. 
filter_canonical_reads_by_mod_base = True 

#Strand to use for training based on alignment to initial masked reference (default: '+-' = no strand filtering)
train_strand = '+-'

#Balance training chunks. May be set to false for testing, otherwise set to true. 
balance_chunks = False

#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base_train = 'S'

#Most similar substituted canonical base you will be comparing against. (default = N, any base)
can_base_train = 'N'

#Range of chunk context to use (in bp) for modified base training (default +/- 0) 
mod_chunk_range = 0

#Chunk context to use (in bp) for modified base training (default +/- 0)
can_chunk_range = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
can_chunk_shift = 0

#working directory 
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240509_moved_functions_test'

#Modified raw data directory 
xna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5'

#Reference fasta in xfasta format for modified dataset
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_10N.fa'
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train.fa'

#Canonical raw data directory 
dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240109_AT_90mer_xr_train_rerun/50f5'

#Reference fasta in xfasta format for canonical dataset
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa'
dna_ref_fasta = xna_ref_fasta

############################################################
#####################Run ##################################

#Step 0 - Initialize
print('Xemora [Status] - Initializing Xemora...')

#Generate directories
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
model_dir = check_make_dir(os.path.join(working_dir,'model'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'modified'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_fastq_dir = check_make_dir(os.path.join(mod_dir,'fastq'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
can_dir = check_make_dir(os.path.join(working_dir,'canonical'))
can_pod_dir = check_make_dir(os.path.join(can_dir,'pod5'))
can_fastq_dir = check_make_dir(os.path.join(can_dir,'fastq'))
can_bam_dir = check_make_dir(os.path.join(can_dir,'bam'))
mod_chunk_meta = check_make_dir(os.path.join(chunk_dir,'mod'))
can_chunk_meta = check_make_dir(os.path.join(chunk_dir,'can'))
misc_dir = check_make_dir(os.path.join(working_dir,'misc_data'))

#Step 0: FASTA to xFASTA conversion
'''
major note here, make sure to add lib/ back later 
'''
if os.path.isfile(os.path.expanduser(xna_ref_fasta)): 
    cmd = 'python xr_fasta2x_rc.py '+os.path.expanduser(xna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora [ERROR] - XNA Reference fasta file not found. Please check file exist or file path.')
    sys.exit()

if os.path.isfile(os.path.expanduser(dna_ref_fasta)): 
    cmd = 'python xr_fasta2x_rc.py '+os.path.expanduser(dna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora [ERROR] - DNA Reference fasta file not found. Please check file exist or file path.')
    sys.exit()
    
#Step 1: Raw input data conversion (pod5, fast5)
print('Xemora [STATUS] - Idenifying raw data input type')
xna_filetype = validate_read_directory(xna_raw_dir)
dna_filetype =validate_read_directory(dna_raw_dir)

print('xna file type', xna_filetype)
mod_merge_cmd = merge_convert_reads(xna_raw_dir, xna_filetype, mod_pod_dir, 'merged')


print('dna file type', dna_filetype)
can_merge_cmd = merge_convert_reads(dna_raw_dir, dna_filetype, can_pod_dir, 'merged')

#Basecalling and initial minimap2 alignment 
#add parameter here to turn this off and let you hard input merged pod5/sam file path

if regenerate_basecall == True: 
    mod_aligned_path = preprocessing_reads(xna_ref_fasta, 'modified', xna_raw_dir, mod_pod_dir, mod_fastq_dir, mod_bam_dir, ref_dir)
    can_aligned_path = preprocessing_reads(dna_ref_fasta, 'canonical', dna_raw_dir, can_pod_dir, can_fastq_dir, can_bam_dir, ref_dir)
else: 
    print('Skipping basecalling')
    mod_pod_dir = working_dir + '/modified/pod5' #path to modified merged pod5 directory, do not include /merged.pod5 
    can_pod_dir = working_dir + '/canonical/pod5' #path to canonical merged pod5 directory, do not include /merged.pod5
    mod_aligned_path = working_dir + '/modified/bam/bam_mm2.sam' #path to modified sam file 
    can_aligned_path = working_dir + '/canonical/bam/bam_mm2.sam' #path to canonical sam file

#Generate dataframe from relevant sam file features 
df_mod = sam_to_df(mod_aligned_path)
df_can = sam_to_df(can_aligned_path)

#Generate tuples of CIGAR operations from CIGAR string 
mod_cigar_tuples = parse_cigar(df_mod['CIGAR'].tolist())
can_cigar_tuples = parse_cigar(df_can['CIGAR'].tolist())


#Ensure XNA labels being used are defined. Extracts list of reference header XNA label and XNA position in the original masked reference
if mod_base_train in ''.join(xna_base_pairs): 
    mod_xna_pos, mod_xna_label = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)), df_mod['Reference Sequence Name'])
    can_xna_pos, can_xna_label = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta)), df_can['Reference Sequence Name'])
else: 
    print('Xemora [ERROR] - Unrecognized XNA found, please check "xna_base_pairs" and "mod_base" in xr_params. Exiting...')
    sys.exit()


#Find xna position on every read.alignment. Output is a list of tuples containing the xna position, readID, strand, and basecalled base at XNA position
print('Xemora [STATUS] - Finding XNA positions for every read')
mod_read_xna_pos_all, mod_no_xna = find_xna_positions(mod_cigar_tuples, df_mod['Sequence'].tolist(), df_mod['Position'].tolist(), mod_xna_pos, mod_xna_label, df_mod['Query Name'].tolist(), df_mod['Flag'].tolist())
can_read_xna_pos_all, can_no_xna = find_xna_positions(can_cigar_tuples, df_can['Sequence'].tolist(), df_can['Position'].tolist(), can_xna_pos, can_xna_label, df_can['Query Name'].tolist(), df_can['Flag'].tolist())


#Filter list of reads by strand
print('Xemora [STATUS] - filtering reads by desired strand')
mod_read_xna_pos = filter_by_strand(mod_read_xna_pos_all, train_strand)
can_read_xna_pos = filter_by_strand(can_read_xna_pos_all, train_strand)

#Filter by modified base label - used for context matching 
print('Xemora [STATUS] - filtering reads by modified base label')
mod_read_xna_pos = filter_by_label(mod_read_xna_pos, mod_base_train)
if filter_canonical_reads_by_mod_base: 
    can_read_xna_pos = filter_by_label(can_read_xna_pos, mod_base_train)

#Data processing summary checkpoint - Generate CSV containing relevant data from the list of tuples (both filtered and unfiltered)
print('Xemora [STATUS] - generating CSV files for XNA position basecall')
xna_pos_csv(mod_read_xna_pos_all, os.path.join(misc_dir,'mod_read_xna_pos_all.csv'))
xna_pos_csv(can_read_xna_pos_all, os.path.join(misc_dir,'can_read_xna_pos_all.csv'))
xna_pos_csv(mod_read_xna_pos, os.path.join(misc_dir,'mod_read_xna_pos.csv'))
xna_pos_csv(can_read_xna_pos, os.path.join(misc_dir,'can_read_xna_pos.csv'))

#Generate a BED file for every read (in filtered list) containing XNA position information 
print('Xemora [STATUS] - generating bed file all reads. Bed file generated post filtering.')
mod_bed_path = generate_bed_file(mod_read_xna_pos, ref_dir,'modified', mod_chunk_range, mod_chunk_shift)
can_bed_path = generate_bed_file(can_read_xna_pos, ref_dir, 'canonical', can_chunk_range, can_chunk_shift)

print('Xemora [STATUS] - Filtering dataset (SAM file) by reads in BAM file. Generating Fasta file')
mod_filter_sam, mod_reads_fasta = sam_filter_fasta_gen(mod_aligned_path, ref_dir, mod_bed_path, 'mod', mod_bam_dir)
can_filter_sam, can_reads_fasta = sam_filter_fasta_gen(can_aligned_path, ref_dir, can_bed_path, 'can', can_bam_dir)

#Realign the reads that fit criteria above to themselves using Minimap2
print('Xemora [STATUS] - Performing 1:1 read alignment using Minimap2')
mod_final_bam = minimap2_realign(mod_filter_sam, mod_reads_fasta)
can_final_bam = minimap2_realign(can_filter_sam, can_reads_fasta)

#Remora stuff
#Step 5: Generate Chunks. 
if regenerate_chunks == True:

    #Chunk generation for for modified chunks
    print('Xemora [STATUS] - Generating chunks for modified base training.')
    mod_chunk_path = os.path.join(chunk_dir,'mod_chunks.npz')
    if os.path.exists(mod_chunk_path):
        try: 
            os.remove(mod_chunk_path) 
        except: 
            print('Xemora [ERROR] - Failed to overwrite existing modified chunk file. Check permissions')
            print('Exiting...')
            sys.exit()

    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,'merged.pod5')+' \
       '+mod_final_bam + ' \
      --output-remora-training-file '+mod_chunk_path+' \
      --focus-reference-positions '+mod_bed_path+' \
      --mod-base '+mod_base_train+' '+mod_base_train+' \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)


    #Chunk generation for canonical chunks
    print('Xemora [STATUS] - Generating chunks for canonical base training.')
    can_chunk_path = os.path.join(chunk_dir,'can_chunks.npz')
    if os.path.exists(can_chunk_path):
        try: 
            os.remove(can_chunk_path) 
        except: 
            print('Xemora [ERROR] - Failed to overwrite existing canonical chunk file. Check permissions')
            print('Exiting...')
            sys.exit()

    cmd = 'remora \
      dataset prepare \
      '+os.path.join(can_pod_dir,'merged.pod5')+' \
       '+can_final_bam + ' \
      --output-remora-training-file '+can_chunk_path+' \
      --focus-reference-positions '+can_bed_path+' \
      --mod-base-control \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)

if remerge_chunks == True: 
    print('Xemora [STATUS] - Merging chunks for training.')

    
    if balance_chunks: 
        cmd = 'remora \
          dataset merge \
          --balance \
          --input-dataset '+mod_chunk_path+' '+chunk_num+'_000 \
          --input-dataset '+can_chunk_path+' '+chunk_num+'_000 \
          --output-dataset '+os.path.join(chunk_dir,'training_chunks.npz')
    else: 
        cmd = 'remora \
          dataset merge \
          --input-dataset '+mod_chunk_path+' '+chunk_num+'_000 \
          --input-dataset '+can_chunk_path+' '+chunk_num+'_000 \
          --output-dataset '+os.path.join(chunk_dir,'training_chunks.npz')
    os.system(cmd)

if gen_model == True:
    print('Xemora [STATUS] - Training model.')
    cmd = 'remora \
      model train \
      '+os.path.join(chunk_dir,'training_chunks.npz')+' \
      --model '+ml_model_path+' \
      --device 0 \
      --output-path '+model_dir+' \
      --overwrite \
      --kmer-context-bases '+kmer_context+' \
      --chunk-context '+chunk_context+' \
      --val-prop '+val_proportion+' \
      --batch-size 100 '# + '\
      #--ext-val ' + ext_dir

    os.system(cmd)

    print('Xemora [Status] - Complete. Saving model to '+model_dir)
"""
need to implement the following: turning the kept reads into a fasta file. 
something to parse through the cigar string to get the xna position in the 
training set. realignment of the data to itself. generation of a merged pod5 and
bam file 
"""
