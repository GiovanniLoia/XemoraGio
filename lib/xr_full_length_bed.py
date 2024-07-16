########################################################################
########################################################################
"""
xr_xfasta2bed.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
########################################################################


from xr_params import *
from xr_tools import *
import sys
import numpy as np
from Bio import SeqIO
import os


input_fasta = sys.argv[1]
output_file = sys.argv[2]
xna_base = sys.argv[3]

def fasta_to_bed(input_fasta, output_file, xna_base):
    with open(output_file, "w") as bed_file:
        for record in SeqIO.parse(os.path.expanduser(input_fasta), "fasta"):
            header = record.id
            seq_length = len(record.seq)
            
            if 'GAP' not in header: 
                x_pos_base = fetch_xna_pos(header)
                x_pos_to_rc =[]
                    
                for x in x_pos_base: 
                     x_base = x[0]


            # Determine the strand
            if x_base == xna_base:
                strand = '+'
                bed_line = f"{header}\t0\t{seq_length}\t{header}\t0\t{strand}\n"
                bed_file.write(bed_line)
            elif x_base == xna_base_rc(xna_base,xna_base_pairs): 
                strand = '-'
                bed_line = f"{header}\t0\t{seq_length}\t{header}\t0\t{strand}\n"
                bed_file.write(bed_line)
            else:
                # If neither xna_base nor its reverse complement are found, you may decide on a default
                print('Xemora [ERROR] - XNA not detected, exiting')
                sys.exit()



fasta_to_bed(input_fasta, output_file, xna_base)
                    
                    

