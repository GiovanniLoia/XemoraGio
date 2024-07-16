########################################################################
########################################################################
"""
xr_merge_chunks.py 

Title: Unpublished work

By: J. Sumabat, J. A. Marchand

Updated: 6/10/23

Convience script to allow the merging of chunks in a directory for model generation
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

############################################################
######################Parameters###########################

#Chunk directory to you want to merge 
chunk_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/P2_chunks_and_models/240617_N44_p2_chunks/N_5000_pod5_41-50/chunks'

#Output directory for chunks 
output_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/P2_chunks_and_models/240617_N44_p2_chunks/N_5000_pod5_41-50'

#Chunk type: can either be 'mod' or 'can' 
chunk_type = 'can'

#Merge method for chunk directory. Options: 'bulk' or 'iterative'
merge_chunk_option = "bulk"

#Verifying output directory is made 
output_dir = check_make_dir(output_dir)

def merge_chunks(chunk_dir, output_dir, chunk_type): 
    """
    merge_chunk will merge a directory of chunks to create a merged chunk for 
    Remora model training
    """
    
    #Calls the list_raw_data method to generate a list of all chunk files in a directory
    chunk_file_names = list_raw_data(chunk_dir)

    #Convert the chunk file names to full file pathways 
    chunk_files = []
    for chunk_name in chunk_file_names:
        chunk_file = os.path.join(chunk_dir, chunk_name)
        chunk_files.append(chunk_file) 
    
    #Start of merge command for all merge types 
    merge_cmd_base = 'remora dataset merge'
    
    # Paths for intermediate and final merged chunk files
    initial_merged_chunk = os.path.join(output_dir, 'initial_merged_chunk.npz')
    temp_merged_chunk_file = os.path.join(output_dir, 'temp_merged_chunk.npz')
    final_merged_chunk = os.path.join(output_dir, chunk_type+'_chunk_merged.npz')
    
    # Validate merge_chunk_option
    if merge_chunk_option not in ["bulk", "iterative"]:
        print("Xemora [ERROR] - Chunk merge method either not specified or incorrect, please verify merge_chunk_option is set to either 'bulk' or 'iterative'")
        sys.exit(1)
        
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
            
            
merge_chunks(chunk_dir, output_dir, chunk_type)
