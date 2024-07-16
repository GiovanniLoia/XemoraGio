########################################################################
"""
xr_basecall_all.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 5/2/23


Note - Reference file need not be the same format for basecalling. 
In fact, you can use an unmasked reference file. 
Goal of a reference file is to just suggest position of modification. 


Change log: 
-Removed calculation of purine/pyrimidine labels. Labels are never used downstream. Not needed. 
    labels also hardcode structural information that is not necessarily relevant 
-Changed variable name "basecall" to "raw_basecall" 
-Changed other variable name "basecall" to "regenerate basecall"
-Fixed pod5 files not overwriting, resulting in odd behavior if you re run analysis
-Fixed chunk files not overwriting, resulting in odd behavior if you re run analysis
-Added print output to terminal 
-Added chunk position ranges
-Added chunk position shifts 
-Added automatic detection of modified base from loading model 
-Moved execution of pod5 merging/fast5 conversion to merge_convert_reads function, rather than merge_pod_command (deprecated)
-Added function to generate fasta output with XNA labels. Note - this function currently depends on sam df intermediate file. Probably bad. 
-Note - Param files has BN, SN as confounding pairs. 
-Added misc file output (sam alignment dataframe now exported)
-Added mapping score filter to xr_params.py. This is a first alignment filter


To do for basecall development: 
-Having figured out what could be wrong with the extraction (raw data), will be good to now check extraction with N base. 

Troubleshooting - Why alignments are poor. 
Modified XNA position looks like its sometimes off by a large amount
Going to try - Higher qscore? 
Why this matters - output reads are going to matter

"""
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


from remora.data_chunks import RemoraRead
from remora.model_util import load_model


from remora.inference import call_read_mods

############################################################
######################Parameters###########################

#Strand to check (+, -, or +-)
check_strand = '+-'

#toggle basecalling 
regenerate_basecall = True

#Print results summary to terminal after generating summary file 
print_results_summary = True

#Write Fasta output file using aligned input file and .tsv per read bc file
write_fasta_output = True

#Range of chunk context to use (in bp) for modified base training (default +/- 0) 
mod_chunk_range = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0

#working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240502_BSn_all_strands_larger_dataset' #the double strand model
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240509_moved_functions_test/bc_ST/S0'

#Raw data directory 
raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/fast5_50-99_basecall' #BSn
#raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240109_AT_90mer_xr_train_rerun/fast5_50-72_basecall' #AT

#Reference fasta in xfasta format
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train.fa'

#10 [0,4] S model
model_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240509_moved_functions_test/model/model_best.pt'


############################################################
#####################Run ##################################

################################################
#Step 0 - Initialize
print('Xemora [STATUS] - Initializing Xemora...')

##Generate directories
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
model_dir = check_make_dir(os.path.join(working_dir,'model'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_fastq_dir = check_make_dir(os.path.join(mod_dir,'fastq'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
mod_chunk_meta = check_make_dir(os.path.join(chunk_dir,'mod'))
misc_dir = check_make_dir(os.path.join(working_dir,'misc_data'))


#Load model and use model 
print('Xemora [STATUS] - Loading Xemora model...')
model, model_metadata = load_model(model_file)
mod_base_infer = model_metadata['mod_bases']
can_base_infer = model_metadata['motif_0']

print('Xemora [STATUS] - Xemora model loaded. XNA base = '+mod_base_infer+' | Context = '+can_base_infer+'.')


#Reference FASTA to xFASTA conversion
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
print('Xemora [STATUS] - Idenifying raw data input type')
xna_filetype = validate_read_directory(raw_dir)

print('Xemora [STATUS] - raw data file type detected: '+ xna_filetype)
mod_merge_cmd = merge_convert_reads(raw_dir, xna_filetype, mod_pod_dir, 'merged')

#List of raw data files in raw directory (meant for batching data in future)
mod_raw_list = list_raw_data(raw_dir)

#Basecalling and initial minimap2 alignment 
if regenerate_basecall == True: 
    mod_aligned_path = preprocessing_reads(mod_raw_list[0], xna_ref_fasta, 'modified', raw_dir, mod_pod_dir, mod_fastq_dir, mod_bam_dir, ref_dir)
else: 
    mod_pod_dir = os.path.join(working_dir, 'preprocess/pod5') #path to modified merged pod5 directory, do not include /merged.pod5 
    mod_aligned_path = os.path.join(working_dir, 'preprocess/bam/bam_mm2.sam') #path to modified sam file 
    print('Xemora [STATUS] - Skipping basecalling. Attempting to use files in '+mod_pod_dir)
    
#Generate dataframe from relevant sam file features 
df_mod = sam_to_df(mod_aligned_path)
df_mod.to_csv(os.path.join(misc_dir,'sam_df.csv'))

#Generate tuples of CIGAR operations from CIGAR string 
mod_cigar_tuples = parse_cigar(df_mod['CIGAR'].tolist())

#Ensure XNA labels being used are defined. 
if mod_base_infer in ''.join(xna_base_pairs): 
    mod_xna_pos, mod_xna_label  = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)), df_mod['Reference Sequence Name'])
else: 
    print('Xemora [ERROR] - Unrecognized XNA found in xFasta reference, please check "xna_base_pairs" and "mod_base" in xr_params')
    sys.exit()

#Find xna position on every read. Output is a list of tuples containing the xna position, readID, strand, and basecalled base at XNA position
print('Xemora [STATUS] - Finding XNA positions for every read')
mod_read_xna_pos_all = find_xna_positions(mod_cigar_tuples, df_mod['Sequence'].tolist(), df_mod['Position'].tolist(), mod_xna_pos, mod_xna_label, df_mod['Query Name'].tolist(), df_mod['Flag'].tolist())

#Filter list of reads by strand
print('Xemora [STATUS] - filtering reads by desired strand')
mod_read_xna_pos = filter_by_strand(mod_read_xna_pos_all, check_strand)

#Filter by modified base label (required for basecalling inference) 
print('Xemora [STATUS] - filtering reads by modified base label')
mod_read_xna_pos = filter_by_label(mod_read_xna_pos, mod_base_infer)

#Generate CSV containing relevant data from the list of tuples (both filtered and unfiltered)
print('Xemora [STATUS] - generating CSV files for XNA position basecall')
xna_pos_csv(mod_read_xna_pos_all, os.path.join(misc_dir,'basecalling_reads_all.csv'))
xna_pos_csv(mod_read_xna_pos, os.path.join(misc_dir,'basecalling_reads_filtered.csv'))

#Generate a BED file for every read (in filtered list) containing XNA position information 
print('Xemora [STATUS] - generating bed file for all the reads. Bed file generated post filtering.')
basecall_bed_path = generate_bed_file(mod_read_xna_pos, ref_dir,'modified', mod_chunk_range, mod_chunk_shift)

print('Xemora [STATUS] - Filtering dataset (SAM file) by reads in BAM file. Generating Fasta file')
mod_filter_sam, mod_reads_fasta = sam_filter_fasta_gen(mod_aligned_path, ref_dir, basecall_bed_path, 'mod', mod_bam_dir) #this is single strand pipeline didnt switch back, need to add conditionals 

#Realign the reads that fit criteria above to themselves using Minimap2
print('Xemora [STATUS] - Performing 1:1 read alignment using Minimap2')
mod_final_bam = minimap2_realign(mod_filter_sam, mod_reads_fasta)

#Remora stuff
#Step 5: Generate Chunks. 
chunk_path = os.path.join(chunk_dir,'basecall_chunks.npz')
if regenerate_chunks == True:
    #Delete existing chunk file in this directory
    
    if os.path.exists(chunk_path):
        try: 
            os.remove(chunk_path) 
        except: 
            print('Xemora [ERROR] - Failed to overwrite existing Chunk file. Check permissions')
            print('Exiting...')
            sys.exit()
    
    print('Xemora [STATUS] - Generating chunks for modified basecalling.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,'merged.pod5')+' \
      '+mod_final_bam+' \
      --output-remora-training-file '+chunk_path+' \
      --focus-reference-positions '+basecall_bed_path+' \
      --mod-base '+mod_base_infer+' '+mod_base_infer+' \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)

#Check that chunk generation worked
if os.path.exists(chunk_path) == False: 
    print('Xemora [ERROR] - No chunk file generated. Check filtering to ensure reads are being used for chunk generation.')
    sys.exit()
    
#Step 6 - Run inference 
try: 
    print('Xemora  [STATUS] - Performing basecalling.')
    cmd = 'remora \
      validate from_remora_dataset \
      '+chunk_path+' \
      --model '+os.path.expanduser(model_file)+' \
      --full-results-filename '+os.path.join(working_dir,'per-read_modifications.tsv')+' \
      --out-file '+os.path.join(working_dir,'summary_modifications.tsv')
    os.system(cmd)

    print('Xemora [STATUS] - Basecalling complete.')
    print('Xemora [STATUS] - Saving per-read modification summaries to: '+os.path.join(working_dir,'per-read_modifications.tsv'))
    print('Xemora [STATUS] - Saving summary result file to: '+os.path.join(working_dir,'summary_modifications.tsv'))
except:
    print('Xemora [ERROR] - Failed to initialize basecalling model. Check logs for details.')
    print('Xemora [ERROR] - Exiting...')
    sys.exit()
    



#Step - Generate additional outputs
#mod_fastq_dir
#This should turn in a general function, but needs rewriting portions
#gen_fasta(per_read_mod_file, bam_input? sam_input?, modified_base_label)
#Currently the raw read info is hard coded. sam_to_df was not really working with other files in the process directory

if write_fasta_output == True: 

    print('Xemora [STATUS] - Generating fasta-x output results file')
    df_mod = sam_to_df(mod_aligned_path)
    df_pr_bc = pd.read_csv(os.path.join(working_dir,'per-read_modifications.tsv'), delimiter= "\t")
    df_pr_bc_mod = df_pr_bc[df_pr_bc['class_pred'] == 1]
    fasta_output_path = os.path.join(working_dir,'xna_fasta.fasta')
    with open(fasta_output_path, "w") as file:
        for i in range(0,len(df_pr_bc_mod)):
            #Get read ID
            rid = df_pr_bc_mod.iloc[i]['read_id']
            
            #Get modified XNA position checked
            pos = df_pr_bc_mod.iloc[i]['read_focus_base']
            
            #Get raw sequence 
            raw_seq = df_mod[df_mod['Query Name']==rid]['Sequence'].iloc[0]
            
            strand_flag = df_mod[df_mod['Query Name']==rid]['Flag'].iloc[0]
            if strand_flag == 16:
                #Reverse aligned sequence
                raw_seq = str(Seq(raw_seq).reverse_complement())
                
                #Reverse index position correction 
                pos=pos-1

            #Modified sequence generation
            mod_seq = raw_seq[:pos] + mod_base_infer + raw_seq[pos+1:]

            #Fasta header for sequence 
            header = '>'+rid+'-XPOS['+mod_base_infer+str(pos)+']'
            file.write(f"{header}\n{mod_seq}\n")
    file.close()
    print('Xemora [STATUS] - Saving fasta-x output result file to: '+fasta_output_path)
    
    
    
#Step 7 
if print_results_summary == True: 
    print('Xemora [STATUS] - Basecall summary:')
    results_summary_path=os.path.join(working_dir,'summary_modifications.tsv')
    # Load the CSV file into a DataFrame
    df = pd.read_csv(results_summary_path,delimiter='\t')

    # Display the headers and their corresponding values from the second row (index 1)
    print(df)
    print('Xemora [STATUS] - Exiting...')


