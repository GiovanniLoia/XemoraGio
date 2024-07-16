########################################################################
########################################################################
"""
xr_mod_chunk_all.py 

Title: Unpublished work

By: J. Sumabat, J. A. Marchand

Updated: 6/10/23

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
from alive_progress import alive_bar
import time
import shutil
from xr_tools  import *
from xr_params import *
from xr_all_context_methods import *

######################Parameters###########################

#Set the amount of POD5/FAST5 files to be merged, inputs can be either 1, any other integer, or 'all'
pod_merge = 100

#Toggle basecalling 
regenerate_basecall = False

#Strand to use for training based on alignment to initial masked reference (default: '+-' = no strand filtering)
train_strand = '+'

#Merge chunk directory into 1 modified or canonical chunk file 
merge_chunk_dir = True

#Merge method for chunk directory. Options: 'bulk' or 'iterative'
merge_chunk_option = "bulk"

#Perform 1:1 read alignment using Mimimap2, if toggled off, assumes path already exists
mm2_realign = True

#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base_train = 'S'

#Most similar substituted canonical base you will be comparing against. (default = N, any base)
can_base_train = 'N'

#Range of chunk context to use (in bp) for modified base training (default +/- 0) 
mod_chunk_range = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0

#working directory 
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/DsPx_chunks_and_models/DsPx_flongle_chunk_subset'

#Modified raw data directory 
xna_raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/240227_BSn_IDT_lib_p2/20240227_1052_P2S-01686-B_PAU96685_48edd116/0-4999pod5' 

#Reference fasta in xfasta format for modified dataset
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/231207_BSn_IDT_Xemora/isoG_isoC.fasta' #BSn fasta

############################################################
#####################Run ##################################

#Step 0 - Initialize
print('Xemora [Status] - Initializing Xemora...')

#Generate directories
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'modified'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_bc_dir = check_make_dir(os.path.join(mod_dir,'basecall'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
mod_chunk_dir = check_make_dir(os.path.join(chunk_dir,'mod'))

def batch_file_processing(pod_file, pod_dir, bc_dir, bam_dir, chunk_dir, xfasta_file_path, datatype): 
    """
    datatype = 'modified or canonical'
    """
    #Create the directory name using the individual pod5 name 
    single_bam_dir = check_make_dir(os.path.join(bam_dir, pod_file))
    if regenerate_basecall:
        print('****************************************************************************************************************************************************')
        #Run the preprocessing function 
        mod_aligned_path = preprocessing_reads_single(pod_file, xfasta_file_path, datatype, pod_dir, bc_dir, single_bam_dir, ref_dir)
    else: 
        print('Skipping basecalling')
        mod_aligned_path = single_bam_dir+'/'+pod_file+'_mm2.sam' #path to modified sam file
    
    #Generate dataframe from relevant sam file features 
    df_mod = sam_to_df(mod_aligned_path)

    #Generate tuples of CIGAR operations from CIGAR string 
    mod_cigar_tuples = parse_cigar(df_mod['CIGAR'].tolist())

    #Ensure XNA labels being used are defined. Extracts list of reference header XNA label and XNA position in the original masked reference
    if mod_base_train in ''.join(xna_base_pairs): 
        mod_xna_pos, mod_xna_label = ref_xnapos_gen(xfasta_file_path, df_mod['Reference Sequence Name'])
    else: 
        print('Xemora [ERROR] - Unrecognized XNA found, please check "xna_base_pairs" and "mod_base" in xr_params. Exiting...')
        sys.exit()

    #Find xna position on every read.alignment. Output is a list of tuples containing the xna position, readID, strand, and basecalled base at XNA position
    print('Xemora [STATUS] - Finding XNA positions for every read')
    mod_read_xna_pos_all, mod_no_xna = find_xna_positions(mod_cigar_tuples, df_mod['Sequence'].tolist(), df_mod['Position'].tolist(), mod_xna_pos, mod_xna_label, df_mod['Query Name'].tolist(), df_mod['Flag'].tolist())

    #Filter list of reads by strand
    print('Xemora [STATUS] - filtering reads by desired strand')
    mod_read_xna_pos = filter_by_strand(mod_read_xna_pos_all, train_strand)
    
    #Filter by modified base label - used for context matching 
    print('Xemora [STATUS] - filtering reads by modified base label')
    mod_read_xna_pos = filter_by_label(mod_read_xna_pos, mod_base_train)
    
    #Data processing summary checkpoint - Generate CSV containing relevant data from the list of tuples (both filtered and unfiltered)
    print('Xemora [STATUS] - generating CSV files for XNA position basecall')
    xna_pos_csv(mod_read_xna_pos_all, os.path.join(single_bam_dir,'mod_read_xna_pos_all.csv'))
    xna_pos_csv(mod_read_xna_pos, os.path.join(single_bam_dir,'mod_read_xna_pos.csv'))
    
    #Generate a BED file for every read (in filtered list) containing XNA position information 
    print('Xemora [STATUS] - generating bed file all reads. Bed file generated post filtering.')
    '''
    major note, i edited this function to add a parameter check which number it batch it did 
    '''
    mod_bed_path = generate_bed_file(mod_read_xna_pos, ref_dir, 'modified', mod_chunk_range, mod_chunk_shift, pod_file)

    print('Xemora [STATUS] - Filtering dataset (SAM file) by reads in BAM file. Generating Fasta file')
    mod_filter_sam, mod_reads_fasta = sam_filter_fasta_gen(mod_aligned_path, ref_dir, mod_bed_path, 'mod', single_bam_dir)

    #Realign the reads that fit criteria above to themselves using Minimap2
    if mm2_realign or not os.path.exists(os.path.join(os.path.dirname(mod_filter_sam), 'final.bam')):
        print('Xemora [STATUS] - Performing 1:1 read alignment using Minimap2')
        mod_final_bam = minimap2_realign(mod_filter_sam, mod_reads_fasta)
    else:
        print('Xemora [STATUS] - Skipping 1:1 read realignment')
        mod_final_bam = os.path.join(os.path.dirname(mod_filter_sam), 'final.bam')
    
    #Step 5: Generate Chunks. 
    if regenerate_chunks == True:

        #Chunk generation for modified chunks
        print('Xemora [STATUS] - Generating chunks for modified base training.')
        mod_chunk_path = os.path.join(chunk_dir,'mod_'+pod_file+'_chunks.npz')
        if os.path.exists(mod_chunk_path):
            try: 
                os.remove(mod_chunk_path) 
            except: 
                print('Xemora [ERROR] - Failed to overwrite existing modified chunk file. Check permissions')
                print('Exiting...')
                sys.exit()

        cmd = 'remora \
          dataset prepare \
          '+os.path.join(pod_dir,pod_file+'.pod5')+' \
           '+mod_final_bam + ' \
          --output-remora-training-file '+mod_chunk_path+' \
          --focus-reference-positions '+mod_bed_path+' \
          --mod-base '+mod_base_train+' '+mod_base_train+' \
          --kmer-context-bases '+kmer_context+' \
          --refine-kmer-level-table '+kmer_table_path+' \
          --refine-rough-rescale '+' \
          --chunk-context '+chunk_context
        os.system(cmd)
        print('Chunk file generated in',mod_chunk_path)

def main(): 
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
        
    #Step 1: Raw input data conversion (pod5, fast5)
    print('Xemora [STATUS] - Identifying raw data input type')
    xna_filetype = validate_read_directory(xna_raw_dir)
    
    #Call process_reads function to either copy or merge pod5 files
    process_reads(xna_raw_dir, xna_filetype, mod_pod_dir, pod_merge)
    
    #List of pod5 files in raw directory (meant for batching data in future)
    mod_raw_list = list_raw_data(mod_pod_dir)
    
    #Starting the remora dataset merge command 
    merge_cmd_base = 'remora dataset merge'
    chunk_files = [] #list to store the chunk names 

    # Paths for intermediate and final merged chunk files
    initial_merged_chunk = os.path.join(chunk_dir, 'initial_merged_chunk.npz')
    temp_merged_chunk_file = os.path.join(chunk_dir, 'temp_merged_chunk.npz')
    final_merged_chunk = os.path.join(chunk_dir, 'mod_chunk_merged.npz')

    # Remove pre-existing merged chunk files
    for file_path in [initial_merged_chunk, temp_merged_chunk_file, final_merged_chunk]:
        if os.path.exists(file_path):
            os.remove(file_path)

    processing_start_time = time.time()
    for pod_file in mod_raw_list:
        pod_name = os.path.splitext(os.path.basename(pod_file))[0]
        
        #Perform batch processing for basecaling, MM2 alignment, bed file generation, etc
        batch_file_processing(pod_name, mod_pod_dir, mod_bc_dir, mod_bam_dir, mod_chunk_dir, os.path.join(ref_dir, 'x' + os.path.basename(xna_ref_fasta)), 'modified')
        
        if merge_chunk_dir == True:
            chunk_file = os.path.join(mod_chunk_dir, 'mod_'+pod_name+'_chunks.npz')
            chunk_files.append(chunk_file)
    processing_finish_time = time.time()

    # Validate merge_chunk_option
    if merge_chunk_option not in ["bulk", "iterative"]:
        print("Xemora [ERROR] - Chunk merge method either not specified or incorrect, please verify merge_chunk_option is set to either 'bulk' or 'iterative'")
        sys.exit(1)

    start_time = time.time()
    if merge_chunk_dir:
        print('Xemora [STATUS] - Merging modifed chunk directory')
        
        # If there is only one chunk file, copy it directly to the final merged chunk file
        if len(chunk_files) == 1:
            shutil.copy(chunk_files[0], final_merged_chunk)
            print('Xemora [STATUS] - Only one chunk file generated. Copied to', final_merged_chunk)

        # Merge all chunk files at once 
        elif merge_chunk_option == "bulk":
            if len(chunk_files) >= 2:
                merge_cmd = merge_cmd_base + ' ' + ' '.join(f'--input-dataset {chunk_file} '+chunk_num+'_000' for chunk_file in chunk_files)
                merge_cmd += ' --output-dataset ' + final_merged_chunk
                os.system(merge_cmd)
                print('Xemora [STATUS] - Merged chunk file generated at', final_merged_chunk)
            else:
                print("Xemora [STATUS] - Not enough chunk files to merge")

        # Merge chunk files iteratively (2 files at a time)
        elif merge_chunk_option == "iterative":
            if len(chunk_files) >= 2:
                # Initial merge of the first two chunk files
                initial_merge_cmd = f"{merge_cmd_base} --input-dataset {chunk_files[0]} "+chunk_num+f"_000 --input-dataset {chunk_files[1]} "+chunk_num+f"_000 --output-dataset {initial_merged_chunk}"
                os.system(initial_merge_cmd)
                print('Initial chunk file made')
                merged_chunk_file = initial_merged_chunk

                # Iteratively merge remaining chunk files
                for chunk_file in chunk_files[2:]:
                    merge_cmd = f"{merge_cmd_base} --input-dataset {merged_chunk_file} "+chunk_num+f"_000 --input-dataset {chunk_file} "+chunk_num+f"_000 --output-dataset {temp_merged_chunk_file} "
                    print(merge_cmd)
                    os.system(merge_cmd)
                    # Move the temp merged chunk to be the new merged chunk
                    merged_chunk_file = temp_merged_chunk_file 

                # Rename the final merged chunk file
                os.rename(merged_chunk_file, final_merged_chunk)
                os.remove(initial_merged_chunk)
                print('Xemora [STATUS] - Merged chunk file generated at', final_merged_chunk)
            else:
                print("Xemora [STATUS] - Not enough chunk files to merge")

    end_time = time.time()
    print('Took', processing_finish_time - processing_start_time, 'to finish processing files')
    print('Took', end_time - start_time, 'to finish merging chunk files')
    
if __name__ == '__main__':
    main()
