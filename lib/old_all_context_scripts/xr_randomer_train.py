########################################################################
########################################################################
"""
xr_train.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23"""
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


############################################################
print('Xemora [Status] - Initializing Xemora...')

#Initialize
#working_dir = os.path.expanduser(sys.argv[1])
#xna_raw_dir = os.path.expanduser(sys.argv[2])
#xna_ref_fasta = os.path.expanduser(sys.argv[3])
#dna_raw_dir = os.path.expanduser(sys.argv[4])
#na_ref_fasta = os.path.expanduser(sys.argv[5])
#ext_dir = "~/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/chunks/training_chunks.npz"

#BSn - Jorge Test
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240427_BSn_only_A_Shift_1'
xna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa'
dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240109_AT_90mer_xr_train_rerun/50f5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa'

#BSn
#working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240428_BSn_only_A_bed_shift_10_jay' #keeping name for convience when model testing 
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240428_BSn_no_filter_all_context'
xna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa'
dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240109_AT_90mer_xr_train_rerun/50f5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_fake_randomer.fa'
'''
#PZ 139 mer
#working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240428_only_G_bed_shift_3_jay'
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240428_PZ_shift_0_no_motif'
xna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/20240215_1810_MN37138_ARS988_4bbd5246/pod5_0-72' #PZ
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/reference/PZ_NB25_xr_Train_3N.fasta' #PZ 
dna_raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/20240216_1817_MN41475_ASE526_f9fc38c7/100_pod5' #PZ
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240216_GC_71merPCR_xr_Train/reference/GC_71mer_xr_Train_3N.fasta' #PZ
'''
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
############################################################
"""
Things to consider: 
Consider changing the basecalling protocol here. If possible, have each pod5 or 
fast5 file be basecalled separately. If looped this way, we can avoid having to 
do bam file merging.

Alternatively, we can run through the entire pipeline as normal. then perform 
alignment using minimap2 on individual fast5 files? Thinking about this a bit more 
this doesn't work. 
"""

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


#Step 1: If input file type is fast5, generate temporary pod5 then basecall, if pod5 start iteratively basecalling 
# 240408: Actually not sure if it should be temporary? Need to send a pod5 file into Remora at the end. Maybe just do a merge at end right before training?

#ls input_folder/*.fast5 | guppy_basecaller --save_path output_folder/basecall --config dna_r9.4.1_450bps.cfg, confirmed this command lets you iterate through, is it even worth converting to pod5???? 

def validate_read_directory(reads_dir):
    directory_exists = os.path.isdir(reads_dir)
    homogenous_files = True
    filetype = ""
    if directory_exists:
        directory_files = os.listdir(reads_dir)
        ext_list = []
        for file in directory_files:
            ext = file.split('.')[-1]
            ext_list.append(ext)
        uniques = list(set(ext_list))
        if (len(uniques) != 1):
            homogenous_files = False
            print('Passed reads directory not homogenous. Filetypes found: {}'.format(uniques))
        else:
            filetype = uniques[0]
            
    return filetype
    
print('Xemora [STATUS] - Idenifying raw data input type')
xna_filetype = validate_read_directory(xna_raw_dir)
dna_filetype =validate_read_directory(dna_raw_dir)

def merge_reads_command(reads_dir, filetype, target_dir_path, file_name):
    # This already assumes the read directory and target directory are valid
    cmd = ""
    subcommand = ""
    os_command = ""
    output_filename = os.path.join(target_dir_path, file_name)+'.pod5'
    
    if filetype == 'fast5':
        subcommand = 'convert fast5'
    elif filetype == 'pod5':
        subcommand = 'merge'
    
    os_command = '{}/*.{}'.format(reads_dir, filetype) #need to account for when users enter their file path as A/B/C/ <-- last slash included
        
    cmd = "pod5 {} {} -o {}".format(subcommand, os_command, output_filename)
    print(cmd)

    return cmd
    
mod_merge_cmd = merge_reads_command(xna_raw_dir, xna_filetype, mod_pod_dir, 'merged')
os.system(mod_merge_cmd)
can_merge_cmd = merge_reads_command(dna_raw_dir, dna_filetype, can_pod_dir, 'merged')
os.system(can_merge_cmd)
    
def list_raw_data(raw_data_dir):
    """
    list_raw_data takes in a file path to a fast5/pod5 dataset and outputs a 
    list containing the file name for every file in the directory.
    
    This will be used to single fast5 file basecalling/processing using Guppy
    """
    file_list = []
    # Iterate over all files in the directory
    for file in os.listdir(raw_data_dir):
            file_list.append(file)
    return file_list
    
    
mod_raw_list = list_raw_data(xna_raw_dir)
can_raw_list = list_raw_data(dna_raw_dir)

def preprocessing_reads(raw_data_file, ref_fasta, dataset_type, raw_dir, pod_dir, fastq_dir, bam_dir): #dataset_type means modified or not, really bad variable name lol 
    """
    preprocessing_read will perform all necessary steps that needs to be batched 
    This will include: initial basecalling, minimap2 alignment, 1:1 bed file 
    generation
    """
    #Step 1: Perform an intiial basecalling using Guppy
    if basecall_pod == True: 
        #cmd = 'ls '+os.path.join(raw_dir, raw_data_file)+' | '+os.path.expanduser(basecaller_path)+' -s '+fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out'
        cmd=os.path.expanduser(basecaller_path)+' -i '+pod_dir+' -s '+fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out ' #-a '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))+' --min_qscore '+str(min_qscore)
        os.system(cmd)
    else:
        print(f'Xemora [STATUS] - Skipping basecalling for {dataset_type} bases')

    #Step 2: Align and transfer info using minimap2 
    # Merge modified bam files
    if os.path.isfile(os.path.join(bam_dir, os.path.basename(bam_dir)) + '.bam') == False or regenerate_bam == True: 

        #Purge directory
        cmd_rmv = 'rm -rf '+bam_dir+'/*'
        os.system(cmd_rmv)
        
        #Set bam directory path 
        bam_path = os.path.join(bam_dir, os.path.basename(bam_dir))
        
        # Merging pass bam files
        cmd_pass = 'samtools merge ' + bam_path + '_pass.bam ' + os.path.join(fastq_dir, 'pass/*.bam -f')
        print(f'Xemora [STATUS] - Merging {dataset_type} PASS BAM files.')
        os.system(cmd_pass)
        
        if merge_fail == True:
            # Merging fail bam files
            print(f'Xemora [STATUS] - Merging {dataset_type} FAIL BAM files.')
            cmd_fail = 'samtools merge ' + bam_path + '_fail.bam ' + os.path.join(fastq_dir, 'fail/*.bam -f')
            os.system(cmd_fail)
            
        #Generate merged bam
        cmd_both = 'samtools merge ' + bam_path + '_all.bam ' + bam_path+'*.bam -f'
        os.system(cmd_both)
        
        #Align using minimap2 and transfer movetable 
        """
        this pipeline might not be correct? one reads are aligned here need to extract xna position on a per read level and generate a bed file 
        
        Need to decide if primary alingment only filtershould happen here
        """
        print('Xemora [STATUS] - Aligning BAM files using minimap2')
        #Flag description: "--score N  0", no deduction for 'N' mismatch; '--secondary no', no secondary aligned reads outputted; '--sam-hit-only', no unaligned reads outputted
        cmd = 'samtools fastq -T "*" '+bam_path+'_all.bam | minimap2 -y -ax map-ont --score-N 0 --secondary no --sam-hit-only --MD '+os.path.join(ref_dir,'x'+os.path.basename(ref_fasta))+ ' - | samtools view -F0x800 -ho ' + bam_path+'_mm2.sam'
        os.system(cmd)
        minimap2_output = bam_path+'_mm2.sam'
    return minimap2_output

mod_aligned_path = preprocessing_reads(mod_raw_list[0], xna_ref_fasta, 'modified', xna_raw_dir, mod_pod_dir, mod_fastq_dir, mod_bam_dir)
can_aligned_path = preprocessing_reads(can_raw_list[0], dna_ref_fasta, 'canonical', dna_raw_dir, can_pod_dir, can_fastq_dir, can_bam_dir)

def sam_to_df(sam_file_path):
    data = []
    # Open the BAM file for reading
    with pysam.AlignmentFile(sam_file_path, "r") as samfile:
        # Iterate over each read in the BAM file
        for read in samfile.fetch():
            # Check if the read is unmapped; skip if true
            if read.is_unmapped:
                continue
            data.append([
            read.query_name,
            read.flag,
            read.reference_name,
            read.reference_start,
            read.cigarstring,
            read.query_sequence
            ])
                
    # Create a DataFrame
    columns = ['Query Name', 'Flag', 'Reference Sequence Name', 'Position', 'CIGAR', 'Sequence']
    df = pd.DataFrame(data, columns=columns)
    return df
    
df_mod = sam_to_df(mod_aligned_path)
df_can = sam_to_df(can_aligned_path)
print(df_can)
                
# CIGAR string parser 
#Parses through the CIGAR string and outputs a tuple, index 0 is the CIGAR operator and index 1 is the amount of bases to use that operator
def parse_cigar(cigar_list):
    pattern = re.compile('([MIDNSHPX=])')
    cigar_tuples_list = []
    for i in range(len(cigar_list)):
        values = pattern.split(cigar_list[i])[:-1]  # Split the string on each code character, and remove the trailing empty string
        parsed = [(values[i+1], int(values[i])) for i in range(0, len(values), 2)]  # pair each operation with its count
        cigar_tuples_list.append(parsed)
    return cigar_tuples_list
    
mod_cigar_tuples = parse_cigar(df_mod['CIGAR'].tolist())
can_cigar_tuples = parse_cigar(df_can['CIGAR'].tolist())

#Generating a list of references sequences and xna positions to be used in the 'hept corrector' and 'corrections' function 
def ref_xnapos_gen(ref_Fasta, rname_list, xna): #edited to use xna containing fasta 
    print('Xemora [STATUS] - Generating reference sequence and xna position lists')
    ref_strings = []
    xna_pos_list = []
    with open(ref_Fasta, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            for j in range(len(rname_list)): 
                if record.id==rname_list[j]: 
                    ref_strings.append(str(record.seq)) #appending the sequence correlated to the reference header 
                    header = r'XPOS\[{}:(\d+)\]'.format(xna) 
                    match = re.search(header, record.id) #looking through the header for 'XPOS\[{}:(\d+)\]'
                    if match:
                        xna_pos_list.append(int(match.group(1))) #extracting and appending the xna from the string above
                    else:
                        print('Xemora [STATUS] - Inputted XNA not detected, please check "mod_base" in xr_params')
    return ref_strings, xna_pos_list
    
def categorize_xnas(xna_base_pairs):
    purines = []
    pyrimidines = []

    for pair in xna_base_pairs:
        # Split the pair into individual bases
        first_base, second_base = pair[0], pair[1]
        
        # Check and add the first base to purines if not already present
        if first_base not in purines:
            purines.append(first_base)
        
        # Check and add the second base to pyrimidines if not already present
        if second_base not in pyrimidines:
            pyrimidines.append(second_base)
    
    return purines, pyrimidines

#Creating lists containing two different XNAs
mod_purines, mod_pyrimidines = categorize_xnas(xna_base_pairs)

if mod_base in mod_purines:
    mod_ref_strings, mod_xna_pos = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)), df_mod['Reference Sequence Name'], mod_base)
    can_ref_strings, can_xna_pos = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta)), df_can['Reference Sequence Name'], mod_base)
elif mod_base in mod_pyrimidines:
    index = mod_pyrimidines.index(mod_base)
    mod_ref_strings, mod_xna_pos = ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)), df_mod['Reference Sequence Name'], mod_purines[index])
    can_ref_strings, can_xna_pos= ref_xnapos_gen(os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta)), df_can['Reference Sequence Name'], mod_purines[index])
else: 
    print('Xemora [ERROR] - Unrecognized XNA found, please check "xna_base_pairs" and "mod_base" in xr_params')
    
def find_xna_positions(cigar_tuples_list, aligned_seq_list, refpos_list, xna_pos_list, readID_list, flags_list):
    """
    Find the positions of modifications (xnas) within multiple read sequences using lists of CIGAR strings and other parameters,
    and return corresponding read IDs and strand information.

    Parameters:
    cigar_tuples_list (list of lists of tuples): List of CIGAR operations for each read.
    aligned_seq_list (list of str): List of sequences of the reads.
    refpos_list (list of int): List of 1-based positions where each read starts aligning in the reference.
    xna_pos_list (list of int): List of 0-based positions of the modifications in the reference sequence.
    readID_list (list of str): List of read IDs corresponding to each read.
    flags_list (list of int): List of SAM flags for determining the strand.

    Returns:
    list of tuples: Each tuple contains the 0-based index of the modification, the read ID, and strand.
    """
    results = []  # To store the results
    for cigar_tuples, aligned_seq, refpos, xna_pos, readID, flag in zip(cigar_tuples_list, aligned_seq_list, refpos_list, xna_pos_list, readID_list, flags_list):
        strand = '-' if (flag & 16) else '+'
        if refpos > xna_pos + 1:
            continue
        current_base = 0
        ref_base = refpos
        #xna_pos_1_based = xna_pos + 1
        if cigar_tuples[0][0] == 'S':
            current_base += cigar_tuples[0][1]
        found = False
        for op, length in cigar_tuples:
            if op == 'S':
                continue
            elif op == 'M' or op == 'I':
                for n in range(length):
                    if ref_base == xna_pos or op == 'I' and ref_base == xna_pos-1: #original if ref_base == xna_pos_1_based or op == 'I' and ref_base == xna_pos_1_based
                        results.append((current_base, readID, strand, aligned_seq[current_base]))
                        found = True
                        break
                    ref_base += 1 if op == 'M' else 0
                    current_base += 1
            elif op == 'D':
                if ref_base <= xna_pos < ref_base + length: #original if ref_base <= xna_pos_1_based < ref_base + length:
                    #results.append((current_base - 1 if current_base > 0 else 0, readID, strand))
                    results.append((None, None, None, None))
                    found = True
                    break
                ref_base += length
            if found:
                break
        if not found:
            results.append((None, None, None, None))

    return results

mod_read_xna_pos_all = find_xna_positions(mod_cigar_tuples, df_mod['Sequence'].tolist(), df_mod['Position'].tolist(), mod_xna_pos, df_mod['Query Name'].tolist(), df_mod['Flag'].tolist())
can_read_xna_pos_all = find_xna_positions(can_cigar_tuples, df_can['Sequence'].tolist(), df_can['Position'].tolist(), can_xna_pos, df_can['Query Name'].tolist(), df_can['Flag'].tolist())

def filter_by_base(xna_positions, specific_base):
    """
    Filters the list of XNA positions to include only those that match a specific base and checks
    the strand based on whether the base is a purine or a pyrimidine.

    Parameters:
    xna_positions (list of tuples): Output from find_xna_positions containing the modification index, read ID, strand, and base.
    specific_base (str): The specific base to filter for (e.g., 'A', 'T', 'G', 'C').

    Returns:
    list of tuples: Filtered list where each tuple contains only entries with the specified base and correct strand orientation.
    """
    purines = ['A', 'G']
    pyrimidines = ['C', 'T']

    filtered_results = []
    for xna_pos, readID, strand, base in xna_positions:
        if base == specific_base:
            if (specific_base in purines and strand == '+') or (specific_base in pyrimidines and strand == '-'):
                filtered_results.append((xna_pos, readID, strand, base))
    return filtered_results
    
mod_read_xna_pos = filter_by_base(mod_read_xna_pos_all, can_base)
can_read_xna_pos = filter_by_base(can_read_xna_pos_all, can_base)

def xna_pos_csv(read_xna_pos, file_path):
    # Count occurrences of A, T, G, C in the entire dataset's fourth column
    all_letters = [row[3] for row in read_xna_pos]
    count = Counter(all_letters)
    
    # Calculate the total of counts to find percentages
    total = sum(count.values()) if count.values() else 1  # Prevent division by zero
    percentages = {key: f"{(value / total) * 100:.2f}%" for key, value in count.items()}
    
    # Ensure each of A, T, G, C has a count and percentage entry, defaulting to 0 and "0.00%"
    for letter in ['A', 'T', 'G', 'C']:
        count.setdefault(letter, 0)
        percentages.setdefault(letter, "0.00%")

    # Prepare the data with additional columns for the first four rows
    updated_data = []
    letters = ['A', 'T', 'G', 'C']  # The order in which we want to display the letters
    
    for i, row in enumerate(read_xna_pos):
        if i < len(letters):
            # Append the letter, its count, and its percentage
            updated_row = row + (letters[i], count[letters[i]], percentages[letters[i]])
        else:
            # Append empty values for rows beyond the first four
            updated_row = row + ('', '', '')
        updated_data.append(updated_row)
    
    # Writing to the CSV file
    with open(file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write each updated tuple to the file
        for row in updated_data:
            writer.writerow(row)

xna_pos_csv(mod_read_xna_pos_all, os.path.join(misc_dir,'mod_read_xna_pos_all.csv'))
xna_pos_csv(can_read_xna_pos_all, os.path.join(misc_dir,'can_read_xna_pos_all.csv'))
xna_pos_csv(mod_read_xna_pos, os.path.join(misc_dir,'mod_read_xna_pos.csv'))
xna_pos_csv(can_read_xna_pos, os.path.join(misc_dir,'can_read_xna_pos.csv'))

def find_pair_name(base, basepairs):
    """ Find which pair the base belongs to and return concatenated pair name. """
    for pair in basepairs:
        if base in pair:
            return ''.join(pair)  # Creates a string like 'PZ' from pair ('P', 'Z')
    return base  # Default to just the base if no pair is found
    
def generate_bed_file(xna_positions, reference_directory, base, basepairs):
    """
    Generates a BED file from the list of XNA positions using a naming convention based on basepairs.

    Parameters:
    xna_positions (list of tuples): Output from find_xna_positions containing positions, read IDs, and strands.
    reference_directory (str): Directory to save the BED file.
    base (str): Base to use for filtering and naming.
    basepairs (list of tuples): List of base pairs for naming convention.
    """
    # Find the appropriate file name using the basepair logic
    file_base_name = find_pair_name(base, basepairs)
    bed_path = os.path.join(reference_directory, f'{file_base_name}.bed')

    purines, pyrimidines = categorize_xnas(basepairs)
    
    with open(bed_path, 'w') as file:
        for xna_pos, readID, strand, basecall in xna_positions:
            if xna_pos is None:
                continue
            chrom = readID
            #range formula is -n, +n+1. gives +/-n base for range
            chromStart = xna_pos#original,no plus one 
            chromEnd = xna_pos + 1 #original is plus one,

            # Determine the correct name to use based on strand and base
            if strand == '+':
                if base in purines:
                    name = base
                elif base in pyrimidines:
                    index = pyrimidines.index(base)
                    name = purines[index]
            elif strand == '-':
                if base in purines:
                    index = purines.index(base)
                    name = pyrimidines[index]
                else:
                    name = base

            score = 0
            file.write(f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t{score}\t{strand}\n")

    return bed_path


mod_full_bed_path = generate_bed_file(mod_read_xna_pos, ref_dir, mod_base, xna_base_pairs)
can_full_bed_path = generate_bed_file(can_read_xna_pos, ref_dir, can_base, standard_base_pairs)

def filter_bed_by_base(bed_path, reference_directory, base):
    """
    Filters a large BED file and outputs a separate file for a specified base.

    Parameters:
    bed_path (str): Path to the input BED file containing mixed base entries.
    output_directory (str): Directory where the output BED file will be saved.
    base (str): The specific base name to filter and output.

    Returns:
    str: Path of the newly created BED file.
    """
    import os

    # Define the output file path for the specified base
    output_bed_path = os.path.join(reference_directory, f'{base}.bed')

    # Open the output file
    with open(output_bed_path, 'w') as outfile:
        # Read the input BED file and write entries that match the specified base
        with open(bed_path, 'r') as infile:
            for line in infile:
                if line.strip():
                    parts = line.split()
                    name = parts[3]  # Assuming 'name' is in the fourth column
                    if name == base:
                        outfile.write(line)

    return output_bed_path
    
mod_bed_path = filter_bed_by_base(mod_full_bed_path, ref_dir, mod_base)
can_bed_path = filter_bed_by_base(can_full_bed_path, ref_dir, can_base)

def sam_filter_fasta_gen(sam_path, ref_dir, bed_path, datatype):
    """
    Filters out reads from a SAM file based on read IDs listed in a BED file, outputs the filtered reads to another SAM file
    and a FASTA file containing sequences from the filtered reads, with each sequence written on a single line.

    Parameters:
    sam_path (str): Path to the input SAM file.
    ref_dir (str): Directory to output the FASTA file.
    bed_path (str): Path to the BED file used to filter reads by read ID.
    datatype (str): Descriptor for the type of data being processed (e.g., 'mod' or 'can') to label the output file.

    Returns:
    tuple: Paths to the output SAM and FASTA files.
    """
    
    output_sam_path = os.path.join(os.path.dirname(sam_path), 'filtered.sam')
    output_fasta_path = os.path.join(ref_dir, f'{datatype}_reads.fasta')
    
    # Extract valid read IDs from the BED file
    valid_read_ids = set()
    with open(bed_path, 'r') as bed_file:
        for line in bed_file:
            if line.strip():
                parts = line.strip().split()
                valid_read_ids.add(parts[0])  # Assuming the first column contains the read IDs

    filtered_records = []  # List to hold SeqRecord objects for FASTA output

    # Filter SAM file and output as SAM
    with pysam.AlignmentFile(sam_path, "r") as sam_file, \
         pysam.AlignmentFile(output_sam_path, "w", template=sam_file) as output_sam:
        for read in sam_file.fetch(until_eof=True):  # Fetch all reads (not relying on indexing)
            if read.query_name in valid_read_ids:
                output_sam.write(read)
                # Create a SeqRecord for each valid read and add to list
                seq_record = SeqRecord(Seq(read.seq), id=read.query_name, description="")
                filtered_records.append(seq_record)

    # Write the FASTA file from the list of SeqRecords, ensuring each sequence is written on a single line
    with open(output_fasta_path, "w") as output_fasta:
        for record in filtered_records:
            output_fasta.write(f">{record.id}\n{str(record.seq)}\n")

    return output_sam_path, output_fasta_path

    
#mod_filter_sam, mod_reads_fasta = sam_filter_fasta_gen(mod_aligned_path, ref_dir, mod_bed_path, 'mod')
#can_filter_sam, can_reads_fasta = sam_filter_fasta_gen(can_aligned_path, ref_dir, can_bed_path, 'can')

def minimap2_realign(sam_file_path, read_fasta):
    output_bam_path = os.path.join(os.path.dirname(sam_file_path), 'final.bam')
    cmd = 'samtools fastq -T "*" '+sam_file_path+' | minimap2 -y -ax map-ont --score-N 0 --secondary no --sam-hit-only --MD '+read_fasta+ ' - | samtools view -F0x800 -bho ' + output_bam_path
    os.system(cmd)
    return output_bam_path

mod_final_bam = minimap2_realign(mod_filter_sam, mod_reads_fasta)
can_final_bam = minimap2_realign(can_filter_sam, can_reads_fasta)


'''
need to add batching to cycle, merging is not cutting it
'''

#Remora stuff
#Step 5: Generate Chunks. 
#os.path.join(xna_raw_dir, mod_raw_list[0]) path to single fast5 file, need to convert 
#os.path.join(dna_raw_dir, can_raw_list[0]), path to single fast5 file, need to convert 
if regenerate_chunks == True:
    
    print('Xemora [STATUS] - Generating chunks for modified base training.')
    '''
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,'merged.pod5')+' \
       '+mod_final_bam + ' \
      --output-remora-training-file '+os.path.join(chunk_dir,'mod_chunks.npz')+' \
      --focus-reference-positions '+mod_bed_path+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    '''
    #no motif allowed
    print('Xemora [STATUS] - Generating chunks for modified base training.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,'merged.pod5')+' \
       '+mod_final_bam + ' \
      --output-remora-training-file '+os.path.join(chunk_dir,'mod_chunks.npz')+' \
      --focus-reference-positions '+mod_bed_path+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context

    os.system(cmd)

#      --focus-reference-positions '+os.path.splitext(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)))[0]+'.bed'+' \
#      '+os.path.join(can_bam_dir,os.path.basename(can_bam_dir))+'.bam'+' \
if regenerate_chunks == True: 
    print('Xemora [STATUS] - Generating chunks for canonical base training.')
    '''
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(can_pod_dir,'merged.pod5')+' \
       '+can_final_bam + ' \
      --output-remora-training-file '+os.path.join(chunk_dir,'can_chunks.npz')+' \
      --focus-reference-positions '+can_bed_path+' \
      --mod-base-control \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)
    '''
    #no motif allowed
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(can_pod_dir,'merged.pod5')+' \
       '+can_final_bam + ' \
      --output-remora-training-file '+os.path.join(chunk_dir,'can_chunks.npz')+' \
      --focus-reference-positions '+can_bed_path+' \
      --mod-base-control \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)

if remerge_chunks == True: 
    print('Xemora [STATUS] - Merging chunks for training.')
    cmd = 'remora \
      dataset merge \
      --balance \
      --input-dataset '+os.path.join(chunk_dir,'mod_chunks.npz')+' '+chunk_num+'_000 \
      --input-dataset '+os.path.join(chunk_dir,'can_chunks.npz')+' '+chunk_num+'_000 \
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
