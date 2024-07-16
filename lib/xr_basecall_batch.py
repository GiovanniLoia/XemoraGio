########################################################################
"""
xr_basecall_batch.py 

Title: Unpublished work

By: J. Sumabat, J. A. Marchand

Updated: 6/13/24

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
import time
import shutil 
from xr_tools  import *
from xr_params import *
from xr_all_context_methods import *


from remora.data_chunks import RemoraRead
from remora.model_util import load_model


from remora.inference import call_read_mods

############################################################
######################Parameters###########################

#Set the amount of POD5/FAST5 files to be merged, inputs can be either 1, any other integer, or 'all'
pod_merge = 'all'

#Toggle pod5 merging
regenerate_pod5 = False #not implemented 

#Toggle basecalling 
regenerate_basecall = True

#Merge chunk directory into 1 modified or canonical chunk file 
merge_chunk_dir = True

#Merge method for chunk directory. Options: 'bulk' or 'iterative' (merge all chunk files at once or one file at a time)
merge_chunk_option = "bulk"

#Perform 1:1 read alignment using Mimimap2, if set to, assumes path already exists
mm2_realign = True

#Strand to check (+, -, or +-)
check_strand = '+-'

#Print results summary to terminal after generating summary file 
print_results_summary = True

#Write Fasta output file using aligned input file and .tsv per read bc file
write_fasta_output = True

#Range of chunk context to use (in bp) for modified base training (default +/- 0) 
mod_chunk_range = 13

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0

#kmer table, variable is overwritten from xr_params import since this is run from lib currrently 
#kmer_table_path = '../models/remora/9mer_10.4.1.csv'

#ml model (ConvLSTM_w_ref.py or Conv_w_ref.py'), variable is overwritten from xr_params import since this is run from lib currrently 
#ml_model_path = '../models/ConvLSTM_w_ref.py'

#working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240502_BSn_all_strands_larger_dataset' #the double strand model
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/model_evaluations/range_tests/SN_sup_+-10'

#Raw data directory 
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/231003_BSn_libv4_FLG114/20231003_1555_MN37138_AQK018_63a7330b/fast5' #BSn xm lib 
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/240227_BSn_IDT_lib_p2/20240227_1052_P2S-01686-B_PAU96685_48edd116/14010-14111pod5' #P2 BSn 
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/240227_BSn_IDT_lib_p2/20240227_1052_P2S-01686-B_PAU96685_48edd116/0-29pod5' #P2 BSn 10 pod5
#raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_AT_90mer_xr_train/20240104_1445_MN41475_ARW614_f45baf6f/fast5'
#raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/300fast5' #BSn single context BSn 
raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5'
#raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240109_AT_90mer_xr_train_rerun/fast5' #AT single se
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/231208_GC_libv4_minion/20231208_1053_MN37138_FAX71315_c0eb1cd1/120fast5'
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/240611_N44_ATGC_lib_p2/20240611_1522_P2S-01686-B_PAU95737_62a852f6/30474-30573pod5' #N44 subset

#Reference fasta
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train.fa'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/231003_BSn_libv4_FLG114/ref_libv2_BS_CxDx-.fa' #BSn xm lib (can be used for GC minion run as well
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/240515_N44_ATGC_lib_flongle/metadata_fasta/ATGC_N44.fa' #N44 fasta
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/231207_BSn_IDT_Xemora/isoG_isoC.fasta' #N44 fasta
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/240227_BSn_IDT_lib_p2/metadata/isoG.fasta'
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Sumabat/240227_BSn_IDT_lib_p2/metadata/isoC.fasta' #P2 Sn

#Xemora model file path 
model_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/SN_balanced_training/model/SN_model.pt' #Balancd SN model 

############################################################
#####################Run ##################################

################################################
#Step 0 - Initialize
print('Xemora [STATUS] - Initializing Xemora...')

##Generate directories
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_bc_dir = check_make_dir(os.path.join(mod_dir,'basecall'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
bc_chunk_dir = check_make_dir(os.path.join(chunk_dir, 'batched_chunks'))

#Load model and use model 
print('Xemora [STATUS] - Loading Xemora model...')
model, model_metadata = load_model(model_file)
mod_base_infer = model_metadata['mod_bases']
can_base_infer = model_metadata['motif_0']

print('Xemora [STATUS] - Xemora model loaded. XNA base = '+mod_base_infer+' | Context = '+can_base_infer+'.')

def batch_file_processing(pod_file, pod_dir, bc_dir, bam_dir, chunk_dir, xfasta_file_path, datatype):
    """
    datatype is always assumed to be modified
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

    #Ensure XNA labels being used are defined. 
    if mod_base_infer in ''.join(xna_base_pairs): 
        mod_xna_pos, mod_xna_label  = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)), df_mod['Reference Sequence Name'])
    else: 
        print('Xemora [ERROR] - Unrecognized XNA found in xFasta reference, please check "xna_base_pairs" and "mod_base" in xr_params')
        sys.exit()

    #Find xna position on every read. Output is a list of tuples containing the xna position, readID, strand, and basecalled base at XNA position
    print('Xemora [STATUS] - Finding XNA positions for every read')
    mod_read_xna_pos_all, mod_no_xna = find_xna_positions(mod_cigar_tuples, df_mod['Sequence'].tolist(), df_mod['Position'].tolist(), mod_xna_pos, mod_xna_label, df_mod['Query Name'].tolist(), df_mod['Flag'].tolist())

    #Filter list of reads by strand
    print('Xemora [STATUS] - filtering reads by desired strand')
    mod_read_xna_pos = filter_by_strand(mod_read_xna_pos_all, check_strand)

    #Filter by modified base label (required for basecalling inference) 
    print('Xemora [STATUS] - filtering reads by modified base label')
    mod_read_xna_pos = filter_by_label(mod_read_xna_pos, mod_base_infer)
    
    #Data processing summary checkpoint - Generate CSV containing relevant data from the list of tuples (both filtered and unfiltered)
    print('Xemora [STATUS] - generating CSV files for XNA position basecall')
    xna_pos_csv(mod_read_xna_pos_all, os.path.join(single_bam_dir,'basecalling_reads_all.csv'))
    xna_pos_csv(mod_read_xna_pos, os.path.join(single_bam_dir,'basecalling_reads_filtered.csv'))
    
    #Generate a BED file for every read (in filtered list) containing XNA position information 
    print('Xemora [STATUS] - generating bed file for all the reads. Bed file generated post filtering.')
    '''
    major note, i edited this function to add a parameter check which number it batch it did 
    '''
    bc_bed_path = generate_bed_file(mod_read_xna_pos, ref_dir,'modified', mod_chunk_range, mod_chunk_shift,pod_file)

    print('Xemora [STATUS] - Filtering dataset (SAM file) by reads in BAM file. Generating Fasta file')
    mod_filter_sam, mod_reads_fasta = sam_filter_fasta_gen(mod_aligned_path, ref_dir, bc_bed_path, 'mod', single_bam_dir) #this is single strand pipeline didnt switch back, need to add conditionals 

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
        bc_chunk_path = os.path.join(chunk_dir,'bc_'+pod_file+'_chunks.npz')
        if os.path.exists(bc_chunk_path):
            try: 
                os.remove(bc_chunk_path) 
            except: 
                print('Xemora [ERROR] - Failed to overwrite existing modified chunk file. Check permissions')
                print('Exiting...')
                sys.exit()

        cmd = 'remora \
          dataset prepare \
          '+os.path.join(pod_dir,pod_file+'.pod5')+' \
           '+mod_final_bam + ' \
          --output-remora-training-file '+bc_chunk_path+' \
          --focus-reference-positions '+bc_bed_path+' \
          --mod-base '+mod_base_infer+' '+mod_base_infer+' \
          --kmer-context-bases '+kmer_context+' \
          --refine-kmer-level-table '+kmer_table_path+' \
          --refine-rough-rescale '+' \
          --chunk-context '+chunk_context
        os.system(cmd)
        print('Chunk file generated in',bc_chunk_path)

    
#Step 6 - Run inference 
def xemora_inference(working_dir, bc_chunk_path, model_file):
    try: 
        print('Xemora  [STATUS] - Performing basecalling.')
        cmd = 'remora \
          validate from_remora_dataset \
          '+bc_chunk_path+' \
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
#This should turn in a general function, but needs rewriting portions
#gen_fasta(per_read_mod_file, bam_input? sam_input?, modified_base_label)
#Currently the raw read info is hard coded. sam_to_df was not really working with other files in the process directory

#Need to add this back sometime, -Jayson
'''
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
'''
def main(): 
    #Step 0: FASTA to xFASTA conversion
    '''
    major note here, make sure to add lib/ back later 
    '''
    if os.path.isfile(os.path.expanduser(xna_ref_fasta)): 
        cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(xna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
        #cmd = 'python xr_fasta2x_rc.py '+os.path.expanduser(xna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
        os.system(cmd)
    else: 
        print('Xemora [ERROR] - XNA Reference fasta file not found. Please check file exist or file path.')
        sys.exit()
        
    #Step 1: Raw input data conversion (pod5, fast5)
    print('Xemora [STATUS] - Identifying raw data input type')
    xna_filetype = validate_read_directory(raw_dir)
    
    #Call process_reads function to either copy or merge pod5 files
    process_reads(raw_dir, xna_filetype, mod_pod_dir, pod_merge)
    
    #List of pod5 files in raw directory (meant for batching data in future)
    mod_raw_list = list_raw_data(mod_pod_dir)
    
    #Starting the remora dataset merge command 
    merge_cmd_base = 'remora dataset merge'
    chunk_files = [] #list to store the chunk names 

    # Paths for intermediate and final merged chunk files
    initial_merged_chunk = os.path.join(chunk_dir, 'initial_merged_chunk.npz')
    temp_merged_chunk_file = os.path.join(chunk_dir, 'temp_merged_chunk.npz')
    final_merged_chunk = os.path.join(chunk_dir, 'bc_chunk_merged.npz')

    # Remove pre-existing merged chunk files
    for file_path in [initial_merged_chunk, temp_merged_chunk_file, final_merged_chunk]:
        if os.path.exists(file_path):
            os.remove(file_path)

    processing_start_time = time.time()
    for pod_file in mod_raw_list:
        pod_name = os.path.splitext(os.path.basename(pod_file))[0]
        batch_file_processing(pod_name, mod_pod_dir, mod_bc_dir, mod_bam_dir, bc_chunk_dir, os.path.join(ref_dir, 'x' + os.path.basename(xna_ref_fasta)), 'modified')
        
        if merge_chunk_dir == True:
            chunk_file = os.path.join(bc_chunk_dir, 'bc_'+pod_name+'_chunks.npz')
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
    
    #Perform Inference
    xemora_inference(working_dir, final_merged_chunk, model_file)

if __name__ == '__main__':
    main()
