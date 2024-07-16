########################################################################
########################################################################
"""
xr_can_chunk_all.py 

Title: Unpublished work

By: J. Sumabat, J. A. Marchand

Updated: 5/29/23

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
filter_canonical_reads_by_mod_base = False

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
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240602_n44_pod5_dir_remora_input_check'

#Canonical raw data directory 
#dna_raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/240515_N27_ATGC_lib_flongle/20240515_1341_MN41475_APP998_377329d7/pod5' #N27 set 
dna_raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/240515_N44_ATGC_lib_flongle/20240515_1335_MN37138_APS278_227b2f57/pod5' #N44 set

#Reference fasta in xfasta format for canonical dataset
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/240515_N27_ATGC_lib_flongle/metadata_fasta/ATGC_N27.fa' #N27 fasta 
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/240515_N44_ATGC_lib_flongle/metadata_fasta/ATGC_N44.fa' #N44 fasta

############################################################
#####################Run ##################################

#Step 0 - Initialize
print('Xemora [Status] - Initializing Xemora...')

#Generate directories
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
model_dir = check_make_dir(os.path.join(working_dir,'model'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
can_dir = check_make_dir(os.path.join(working_dir,'canonical'))
can_pod_dir = check_make_dir(os.path.join(can_dir,'pod5'))
can_fastq_dir = check_make_dir(os.path.join(can_dir,'fastq'))
can_bam_dir = check_make_dir(os.path.join(can_dir,'bam'))
can_chunk_meta = check_make_dir(os.path.join(chunk_dir,'can'))
misc_dir = check_make_dir(os.path.join(working_dir,'misc_data'))

#Step 0: FASTA to xFASTA conversion
'''
major note here, make sure to add lib/ back later 
'''
if os.path.isfile(os.path.expanduser(dna_ref_fasta)): 
    cmd = 'python xr_fasta2x_rc.py '+os.path.expanduser(dna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora [ERROR] - DNA Reference fasta file not found. Please check file exist or file path.')
    sys.exit()
    
#Step 1: Raw input data conversion (pod5, fast5)
print('Xemora [STATUS] - Idenifying raw data input type')
dna_filetype =validate_read_directory(dna_raw_dir)

print('dna file type', dna_filetype)
can_merge_cmd = merge_convert_reads(dna_raw_dir, dna_filetype, can_pod_dir, 'merged')

#Basecalling and initial minimap2 alignment 
#add parameter here to turn this off and let you hard input merged pod5/sam file path

if regenerate_basecall == True: 
    can_aligned_path = preprocessing_reads(dna_ref_fasta, 'canonical', dna_raw_dir, can_pod_dir, can_fastq_dir, can_bam_dir, ref_dir)
else: 
    print('Skipping basecalling')
    can_pod_dir = working_dir + '/canonical/pod5' #path to canonical merged pod5 directory, do not include /merged.pod5
    can_aligned_path = working_dir + '/canonical/bam/bam_mm2.sam' #path to canonical sam file

#Generate dataframe from relevant sam file features 
df_can = sam_to_df(can_aligned_path)

#Generate tuples of CIGAR operations from CIGAR string 
can_cigar_tuples = parse_cigar(df_can['CIGAR'].tolist())


#Ensure XNA labels being used are defined. Extracts list of reference header XNA label and XNA position in the original masked reference
if mod_base_train in ''.join(xna_base_pairs): 
    can_xna_pos, can_xna_label = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta)), df_can['Reference Sequence Name'])
else: 
    print('Xemora [ERROR] - Unrecognized XNA found, please check "xna_base_pairs" and "mod_base" in xr_params. Exiting...')
    sys.exit()

#Find xna position on every read.alignment. Output is a list of tuples containing the xna position, readID, strand, and basecalled base at XNA position
print('Xemora [STATUS] - Finding XNA positions for every read')
can_read_xna_pos_all, can_no_xna = find_xna_positions(can_cigar_tuples, df_can['Sequence'].tolist(), df_can['Position'].tolist(), can_xna_pos, can_xna_label, df_can['Query Name'].tolist(), df_can['Flag'].tolist())

#Filter list of reads by strand
print('Xemora [STATUS] - filtering reads by desired strand')
can_read_xna_pos = filter_by_strand(can_read_xna_pos_all, train_strand)

#Filter by modified base label - used for context matching 
print('Xemora [STATUS] - filtering reads by modified base label')
if filter_canonical_reads_by_mod_base: 
    can_read_xna_pos = filter_by_label(can_read_xna_pos, mod_base_train)

#Data processing summary checkpoint - Generate CSV containing relevant data from the list of tuples (both filtered and unfiltered)
print('Xemora [STATUS] - generating CSV files for XNA position basecall')
xna_pos_csv(can_read_xna_pos_all, os.path.join(misc_dir,'can_read_xna_pos_all.csv'))
xna_pos_csv(can_read_xna_pos, os.path.join(misc_dir,'can_read_xna_pos.csv'))

#Generate a BED file for every read (in filtered list) containing XNA position information 
print('Xemora [STATUS] - generating bed file all reads. Bed file generated post filtering.')
can_bed_path = generate_bed_file(can_read_xna_pos, ref_dir, 'canonical', can_chunk_range, can_chunk_shift)

print('Xemora [STATUS] - Filtering dataset (SAM file) by reads in BAM file. Generating Fasta file')
can_filter_sam, can_reads_fasta = sam_filter_fasta_gen(can_aligned_path, ref_dir, can_bed_path, 'can', can_bam_dir)

#Realign the reads that fit criteria above to themselves using Minimap2
print('Xemora [STATUS] - Performing 1:1 read alignment using Minimap2')
can_final_bam = minimap2_realign(can_filter_sam, can_reads_fasta)

#Remora stuff
#Step 5: Generate Chunks. 
if regenerate_chunks == True:

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


    print('Canonical chunks generated in',can_chunk_path)
"""
need to implement the following: turning the kept reads into a fasta file. 
something to parse through the cigar string to get the xna position in the 
training set. realignment of the data to itself. generation of a merged pod5 and
bam file 
"""
