########################################################################
"""
xr_all_context_functions.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 5/29/24

xr_all_context_methods.py contains the shared functions for xr_train_all.py &
xr_basecall_all.py. This allows for changes in the functions to apply to both 
scripts (do not need to copy paste changes). 


Note - Reference file need not be the same format for basecalling. 
In fact, you can use an unmasked reference file. 
Goal of a reference file is to just suggest position of modification. 


Change log: 
- 5/29: find_xna_positions function was changed to remove reads where XNA position wasn't found. Number of reads removed is now a print statement 
- 5/7: in 'generate_bed_file', added chunk range and shifts as inputs, instead of calling global variables 
- 5/7: added reference directory as an input into function 
-Removed calculation of purine/pyrimidine labels. Labels are never used downstream. Not needed. 
    labels also hardcode structural information that is not necessarily relevant 
-Changed variable name "basecall" to "raw_basecall" 
-Changed other variable name "basecall" to "regenerate basecall"
-Fixed pod5 files not overwriting, resulting in odd behavior if you re run analysis
-Fixed chunk files not overwriting, resulting in odd behavior if you re run analysis
-Added print output to terminal 
-Added chunk position ranges
-Added chunk position shifts 
-Moved execution of pod5 merging/fast5 conversion to merge_convert_reads function, rather than merge_pod_command (deprecated)
-Added versbose output flag 
"""
########################################################################
########################Imports####################################

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
import shutil
from xr_tools  import *
from xr_params import *

############################################################
#####################Functions##############################

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
            print('Xemora [STATUS] - Passed reads directory not homogenous. Filetypes found: {}'.format(uniques))
        else:
            filetype = uniques[0]
            
    return filetype
    
def merge_convert_reads(reads_dir, filetype, target_dir_path, file_index, files_to_process):
    """
    merge_convert_reads takes in a fast5/pod5 directory, the filetype in that directory,
    an output directory, and desired output name to generate and run a command to 
    generate a single merged pod5 file in the desired target directory.
    """
    # Initialize subcommand and os_command variables
    subcommand = ""
    os_command = ""
    # Create the output filename using the provided file index
    output_filename = os.path.join(target_dir_path, f"{file_index}.pod5")
    
    # Determine the subcommand based on the file type
    if filetype == 'fast5':
        subcommand = 'convert fast5'
    elif filetype == 'pod5':
        subcommand = 'merge'
    
    # Combine the list of files to process into a single command string
    os_command = ' '.join(files_to_process)

    # Create the final command for merging/converting files, always using force-overwrite
    cmd = f"pod5 {subcommand} {os_command} -o {output_filename} --force-overwrite"
    # Execute the command
    os.system(cmd)

    return True
    
def single_convert_copy_reads(reads_dir, filetype, target_dir_path):
    """
    single_convert_copy takes in a fast5/pod5 directory and does either single
    fast5 to pod5 conversion or copies the individual pod5 files to the target_dir_path
    directory. The output files are renamed sequentially.
    """
    # Initialize file counter for sequential naming
    file_counter = 1
    # Loop through all files in the reads directory
    for filename in os.listdir(reads_dir):
        source_file_path = os.path.join(reads_dir, filename)  # Full path to the source file
        
        # Check if the file is a fast5 file
        if filetype == 'fast5' and filename.endswith('.fast5'):
            # Create the target file path with sequential naming
            target_file_path = os.path.join(target_dir_path, f"{file_counter}.pod5")
            # Create the command to convert the fast5 file to a pod5 file
            cmd = f"pod5 convert fast5 {source_file_path} -o {target_file_path} --force-overwrite"
            # Execute the command
            os.system(cmd)
 
            # Increment the file counter
            file_counter += 1
        # Check if the file is a pod5 file
        elif filetype == 'pod5' and filename.endswith('.pod5'):
            # Create the target file path with sequential naming
            target_file_path = os.path.join(target_dir_path, f"{file_counter}.pod5")
            # Copy the pod5 file to the target directory
            shutil.copy(source_file_path, target_file_path)
            # Increment the file counter
            file_counter += 1

def process_reads(reads_dir, filetype, target_dir_path, pod_merge):
    """
    process_reads manages the conversion and merging process based on the pod_merge variable.
    """
    # Handle the case where pod_merge is 1 (single convert/copy)
    if pod_merge == 1:
        print('Xemora [STATUS] - Copying/converting files')
        single_convert_copy_reads(reads_dir, filetype, target_dir_path)
        print('Xemora [STATUS] - Files copied/converted')
    # Handle the case where pod_merge is 'all' (merge all files)
    elif pod_merge == 'all':
        print('Xemora [STATUS] - Merging all files')
        merge_convert_reads(reads_dir, filetype, target_dir_path, 1, 
                            [os.path.join(reads_dir, f) for f in os.listdir(reads_dir) if f.endswith(f'.{filetype}')])
    # Handle the case where pod_merge is an integer (batch processing)
    else:
        print('Xemora [STATUS] - Merging', pod_merge, 'files per batch')
        try:
            num_files = int(pod_merge)  # Convert pod_merge to an integer
            files = [os.path.join(reads_dir, f) for f in os.listdir(reads_dir) if f.endswith(f'.{filetype}')]  # List of files to process
            total_files = len(files)  # Total number of files
            # Process files in batches
            for i in range(0, total_files, num_files):
                batch_files = files[i:i+num_files]  # Get the current batch of files
                file_index = i // num_files + 1  # Calculate the file index for naming
                merge_convert_reads(reads_dir, filetype, target_dir_path, file_index, batch_files)  # Call the merge/convert function
        except ValueError:
            # Print an error message if pod_merge is not a valid integer
            print("Invalid value for pod_merge. It must be 1, 'all', or an integer.")

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

def preprocessing_reads(ref_fasta, dataset_type, raw_dir, pod_dir, fastq_dir, bam_dir, ref_dir): #dataset_type means modified or not, really bad variable name lol 
    """
    preprocessing_read will perform all necessary steps that needs to be batched 
    This will include: initial basecalling, minimap2 alignment, 1:1 bed file 
    generation
    
    add parameters to check whether dorado or guppy 
    """
    #Step 1: Perform an intiial basecalling using Guppy
    if basecall_pod == True: 
        cmd=os.path.expanduser(basecaller_path)+' -i '+pod_dir+' -s '+fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out '
        os.system(cmd)
    else:
        print(f'Xemora [STATUS] - Skipping basecalling for {dataset_type} bases')
    '''
    #Step 1: Perform an intiial basecalling using Dorado
    if basecall_pod == True: 
        if dorado_model == '':
            cmd = "{} basecaller {} --no-trim --emit-moves --min-qscore {} {} > {}".format(dorado_path, 'sup', min_qscore, os.path.join(raw_dir, 'bc.pod5'), os.path.join(bc_dir, raw_data_file+'.bam'))
            os.system(cmd)
        else: 
            cmd = "{} basecaller {} --no-trim --emit-moves --min-qscore {} {} > {}".format(dorado_path, dorado_model, min_qscore, os.path.join(raw_dir, raw_data_file+'.pod5'), os.path.join(bc_dir, raw_data_file+'.bam'))
            os.system(cmd)
    else:
        print(f'Xemora [STATUS] - Skipping basecalling for {dataset_type} bases')
    '''
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
            print(f'Xemora [STATUS] - Param merge_fail set to True. Merging {dataset_type} FAIL BAM files.')
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
    
def preprocessing_reads_single(raw_data_file, ref_fasta, dataset_type, raw_dir, bc_dir, bam_dir, ref_dir): #dataset_type means modified or not, really bad variable name lol 
    """
    preprocessing_read will perform all necessary steps that needs to be batched 
    This will include: initial basecalling, minimap2 alignment, 1:1 bed file 
    generation
    """
    #Step 1: Perform an intiial basecalling using Dorado
    if basecall_pod == True: 
        if dorado_model == '':
            cmd = "{} basecaller {} --no-trim --emit-moves --min-qscore {} {} > {}".format(dorado_path, 'sup', min_qscore, os.path.join(raw_dir, raw_data_file+'.pod5'), os.path.join(bc_dir, raw_data_file+'.bam'))
            os.system(cmd)
        else: 
            cmd = "{} basecaller {} --no-trim --emit-moves --min-qscore {} {} > {}".format(dorado_path, dorado_model, min_qscore, os.path.join(raw_dir, raw_data_file+'.pod5'), os.path.join(bc_dir, raw_data_file+'.bam'))
            os.system(cmd)
    else:
        print(f'Xemora [STATUS] - Skipping basecalling for {dataset_type} bases')
        
    #Align using minimap2 and transfer movetable 
    print('Xemora [STATUS] - Aligning BAM files using minimap2')
    
    #Flag description: "--score N  0", no deduction for 'N' mismatch; '--secondary no', no secondary aligned reads outputted; '--sam-hit-only', no unaligned reads outputted
    cmd = 'samtools fastq -T "*" '+os.path.join(bc_dir, raw_data_file+'.bam')+ ' | minimap2 -y -ax map-ont --score-N 0 --secondary no --sam-hit-only --MD '+ref_fasta+ ' - | samtools view -F0x800 -ho ' + bam_dir+'/'+raw_data_file+'_mm2.sam'
    print(cmd)
    os.system(cmd)
    minimap2_output = bam_dir+'/'+raw_data_file+'_mm2.sam'
    return minimap2_output
    
def sam_to_df(sam_file_path):
    data = []
    # Open the SAM file for reading
    with pysam.AlignmentFile(sam_file_path, "r") as samfile:
        # Iterate over each read in the BAM file
        for read in samfile.fetch():
            # Check if the read is unmapped; skip if true
            if read.is_unmapped:
                continue
            if read.mapping_quality < min_map_score: 
                continue
            data.append([
            read.query_name,
            read.flag,
            read.reference_name,
            read.reference_start,
            read.cigarstring,
            read.query_sequence,
            read.mapping_quality
            ])
                
    # Create a DataFrame
    columns = ['Query Name', 'Flag', 'Reference Sequence Name', 'Position', 'CIGAR', 'Sequence','Mapping Quality']
    df = pd.DataFrame(data, columns=columns)
    return df
                
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
    

#Generating a list of references sequences and xna positions to be used in the 'hept corrector' and 'corrections' function 
def ref_xnapos_gen(ref_fasta, rname_list): #edited to use xna containing fasta 
    print('Xemora [STATUS] - Extracting XNA position labels from xfasta headers for mapping XNA positions to reference sequence')
    pos_list = []
    label_list = []
    
    #Open reference fasta to get a list of possible masked references
    try: 
        with open(ref_fasta, "r") as rfasta:
            #For each possible reference
            for record in SeqIO.parse(rfasta, "fasta"):
                #Check if this reference had aligned to one of the reads from first-pass alignment
                for j in range(len(rname_list)): 
                    #If reference header name was found in an alignment sequence
                    if record.id==rname_list[j]: 
                        #Hardcoded assumption - 1 XNA per header
                        #Extract XNA position from header and xna label
                        xna_pos = fetch_xna_pos(rname_list[j])[0]
                        pos_list.append(int(xna_pos[1]))
                        label_list.append(str(xna_pos[0]))
    except:
        print('Xemora [ERROR] - Improperly formatted xFasta file for training. Failed to extract XNA labels from reference header. Exiting...')
        sys.exit()

    return pos_list, label_list
    

def find_xna_positions(cigar_tuples_list, aligned_seq_list, refpos_list, xna_pos_list, xna_label_list, readID_list, flags_list):
    """
    Find the positions of modifications (xnas) within multiple read sequences using lists of CIGAR strings and other parameters,
    and return corresponding read IDs and strand information. 

    Parameters:
    cigar_tuples_list (list of lists of tuples): List of CIGAR operations for each read.
    aligned_seq_list (list of str): List of sequences of the reads.
    refpos_list (list of int): List of 1-based positions where each read starts aligning in the reference.
    xna_pos_list (list of int): List of 0-based positions of the modifications in the reference sequence.
    xna_label_list (list of string): List of the modification label, as determined from the xFasta header at the given position. 
    readID_list (list of str): List of read IDs corresponding to each read.
    flags_list (list of int): List of SAM flags for determining the strand.

    Returns:
    list of tuples: Each tuple contains the 0-based index of the modification, the read ID, strand, and strand corrected XNA label. Note that XNA label returned is the strand alignment corrected label. 
    """
    
    results = []
    no_xna = 0
    for cigar_tuples, aligned_seq, ref_pos, xna_pos, xna_label, readID, flag in zip(cigar_tuples_list, aligned_seq_list, refpos_list, xna_pos_list, xna_label_list, readID_list, flags_list):
        strand = '-' if (flag & 16) else '+'
        if ref_pos > xna_pos + 1:
            continue
        current_base = 0
        xna_label_rc = xna_base_rc(xna_label, xna_base_pairs)
        
        if cigar_tuples[0][0] == 'S':
            current_base += cigar_tuples[0][1]
        found = False
        for op, length in cigar_tuples:
            if op == 'S':
                continue
            elif op == 'M' or op == 'I':
                for n in range(length):
                    if ref_pos == xna_pos or (op == 'I' and ref_pos == xna_pos - 1):
                        if strand == '-':
                            aligned_seq_current_base_rc = xna_base_rc(aligned_seq[current_base], xna_base_pairs + standard_base_pairs)
                            results.append((len(aligned_seq)-current_base, readID, strand, aligned_seq_current_base_rc,xna_label_rc ))
                        else:
                            results.append((current_base, readID, strand, aligned_seq[current_base],xna_label))
                        found = True
                        break
                    ref_pos += 1 if op == 'M' else 0
                    current_base += 1
            elif op == 'D':
                if ref_pos <= xna_pos < ref_pos + length:
                    if strand == '-':
                        aligned_seq_current_base_rc = xna_base_rc(aligned_seq[current_base-1], xna_base_pairs + standard_base_pairs)
                        results.append((len(aligned_seq)-current_base-1 if current_base > 0 else 0, readID, strand, aligned_seq_current_base_rc,xna_label_rc ))
                    else:
                        results.append((current_base - 1 if current_base > 0 else 0, readID, strand, aligned_seq[current_base-1],xna_label))
                    found = True
                    break
                ref_pos += length
            if found:
                break

        if not found:
            no_xna += 1
    print(no_xna, 'reads removed, XNA position was not found')
    return results, no_xna

def filter_by_strand(xna_positions, desired_strand):
    """
    Filters the list of XNA positions to only include reads that match the desired
    strand 
    """
    filtered_results = []
    for xna_pos, readID, strand, base, xna_label in xna_positions:
        if strand in desired_strand: 
            filtered_results.append((xna_pos,readID, strand, base,xna_label))
    return filtered_results
    
def filter_by_label(xna_positions, mod_base_label):
    """
    Filters the list of XNA positions to only include reads that have modified base label
    """
    filtered_results = []
    for xna_pos, readID, strand, base, xna_label in xna_positions:
        #add a parameter here in xr_params that checks for 'both' for more complex applications, for now, making a hard coded variable to check for one strand
        if xna_label == mod_base_label:
            filtered_results.append((xna_pos,readID, strand, base,xna_label))
    return filtered_results

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
            

def generate_bed_file(xna_position_info, reference_directory, chunk_label, chunk_range, chunk_shift, batch_number):
    """
    Generates a BED file from the list of XNA positions using a naming convention based on basepairs.

    Parameters:
    xna_positions_info (list of tuples): Output from find_xna_positions containing positions, read IDs, strands, and xna_label
    reference_directory (str): Directory to save the BED file.
    base (str): Base to use for filtering and naming.
    basepairs (list of tuples): List of base pairs for naming convention.
    Note - xna base here just refers to the location of training on a given read. Canonical reads will use this label too. 
    chunk_label = 'canonical' or 'modified'
    batch_number = set of data used to generate the bed file
    """
    # Find the appropriate file name using the basepair logic
    if chunk_label =='modified':
        file_base_name = 'XY'
        #chunk_range = mod_chunk_range
        #chunk_shift = mod_chunk_shift
    if chunk_label =='canonical':
        file_base_name = 'ATGC'
        #chunk_range = can_chunk_range
        #chunk_shift = can_chunk_shift

    #Automatically generate path for output bed file
    bed_path = os.path.join(reference_directory, f'{file_base_name}_{batch_number}.bed')

    #Write bed file line by line
    with open(bed_path, 'w') as file:
        for xna_pos, readID, strand, raw_basecall, xna_base_label in xna_position_info:
            if xna_pos is None:
                continue

            #range formula is -n, +n+1. gives +/-n base for range
            
            #Start position of base of interest
            chromStart = xna_pos - chunk_range + chunk_shift#original,no plus one 
            
            #End position of base of interest
            chromEnd = xna_pos + chunk_range + chunk_shift + 1 #original is plus one,

            #Label to use in bed file 
            if chunk_label == 'modified':
                base_label = xna_base_label
            if chunk_label == 'canonical': 
                base_label = raw_basecall

            #Required bed file column
            score = 0
            
            #For double-alignment strategy, bed file strand will always be (+) 
            strand_p = '+'
            
            #file.write(f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t{score}\t{strand}\n")
            file.write(f"{readID}\t{chromStart}\t{chromEnd}\t{base_label}\t{score}\t{strand_p}\n")

    return bed_path

def sam_filter_fasta_gen(sam_path, ref_dir, bed_path, datatype, bam_dir):
    """
    added some RC stuff here
    --Jayson
    
    
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
    output_sam_path = os.path.join(bam_dir, 'filtered.sam')
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
        for read in sam_file.fetch(until_eof=True):
            if read.query_name in valid_read_ids:
                # Check the flag for reverse complement
                if read.flag & 16:  # Check if the 0x10 bit is set, indicating reverse strand
                    seq = Seq(read.seq).reverse_complement()
                else:
                    seq = Seq(read.seq)
                
                output_sam.write(read)
                # Create a SeqRecord for each valid read and add to list
                seq_record = SeqRecord(seq, id=read.query_name, description="")
                filtered_records.append(seq_record)

    # Write the FASTA file from the list of SeqRecords, ensuring each sequence is written on a single line
    with open(output_fasta_path, "w") as output_fasta:
        for record in filtered_records:
            output_fasta.write(f">{record.id}\n{str(record.seq)}\n")

    return output_sam_path, output_fasta_path

def minimap2_realign(sam_file_path, read_fasta):
    output_bam_path = os.path.join(os.path.dirname(sam_file_path), 'final.bam')
    cmd = 'samtools fastq -T "*" '+sam_file_path+' | minimap2 -y -ax map-ont --score-N 0 --secondary no --sam-hit-only --MD '+read_fasta+ ' - | samtools view -F0x800 -bho ' + output_bam_path
    print(cmd)
    print('******************************************')
    os.system(cmd)
    print('******************************************')
    #os.system(cmd)
    return output_bam_path

