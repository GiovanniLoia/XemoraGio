import os
import glob
import sys
from pathlib import Path
from xr_tools  import *
from xr_params import *


############################################################
print('Xemora [Status] - Initializing Xemora...')

#Initialize
working_dir = os.path.expanduser(sys.argv[1])
xna_raw_dir = os.path.expanduser(sys.argv[2])
dna_raw_dir = os.path.expanduser(sys.argv[3])
xna_ref_fasta = os.path.expanduser(sys.argv[4])
dna_ref_fasta = os.path.expanduser(sys.argv[5])


#Assigning Variables to directories that should have been generated during train
ref_dir = os.path.join(working_dir,'references')

#Modified Directories
mod_dir = os.path.join(working_dir,'modified')
mod_pod_dir = os.path.join(mod_dir,'pod5')
mod_bam_dir = os.path.join(mod_dir,'bam')

#Canonical Directories
can_dir = os.path.join(working_dir,'canonical')
can_pod_dir = os.path.join(can_dir,'pod5')
can_bam_dir = os.path.join(can_dir,'bam')

mod_bam_path = os.path.join(mod_bam_dir, os.path.basename(mod_bam_dir))
can_bam_path = os.path.join(can_bam_dir, os.path.basename(can_bam_dir))

#Check to see if required files exist 

#pod5 checks 
if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5')==False:
    print('Xemora [ERROR] - Modified merged pod5 file does not exist. Did you run xr_train?')
    sys.exit()
    
if os.path.isfile(os.path.join(can_pod_dir,os.path.basename(dna_raw_dir))+'.pod5')==False:
    print('Xemora [ERROR] - Canonical merged pod5 file does not exist. Did you run xr_train?')
    sys.exit()
    
#Minimap2 checks
if os.path.isfile(os.path.join(mod_bam_dir,'bam_sorted.bam')) == False: 
    print('Xemora [ERROR] - Modified BAM file from Minimap2 does not exist. Did you run xr_train?')
    sys.exit()
    
if os.path.isfile(os.path.join(can_bam_dir,'bam_sorted.bam')) == False: 
    print('Xemora [ERROR] - Canonical BAM file from Minimap2 does not exist. Did you run xr_train?')
    sys.exit()
    
#Bed file checks 
if os.path.isfile(os.path.join(ref_dir,mod_base+'.bed')) == False:
    print('Xemora [ERROR] - Modified bed file does not exist. Did you run xr_train?')
    sys.exit()

if os.path.isfile(os.path.join(ref_dir,can_base+'.bed')) == False:
    print('Xemora [ERROR] - Canonical bed file does not exist. Did you run xr_train?')
    sys.exit()

#Indexing bam files 

#Remora Raw Signal Command

cmd = 'remora \
  analyze plot ref_region \
  --pod5-and-bam '+ os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5'+' '+os.path.join(mod_bam_dir,'bam_sorted.bam')+' \
  --pod5-and-bam '+ os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5'+' '+os.path.join(can_bam_dir,'bam_sorted.bam')+'\
  --ref-regions '+os.path.join(ref_dir,mod_base+'_full.bed')+' \
  --highlight-ranges '+os.path.join(ref_dir,mod_base+'.bed')+   ' \
  --refine-kmer-level-table '+kmer_table_path+' \
  --refine-rough-rescale \
  --log-filename '+ os.path.join(working_dir, 'log.txt')
  
os.system(cmd)
