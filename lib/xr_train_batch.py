########################################################################
########################################################################
"""
xr_train_batch.py 

Title: Unpublished work

By: J. Sumabat, J. A. Marchand

Updated: 6/10/23

To do 6/25/24: need to have this script get called with parameters that can be 
set to generate chunks if needed or can have the chunks generated prior as inputs

6/27: to do - make sure to add a file path verification check
"""

########################################################################
########################################################################

import os
import glob
import sys
from pathlib import Path
import torch
from xr_tools  import *
from xr_params import *
from xr_all_context_methods import *

############################################################
######################Parameters###########################

#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base_train = 'B'

#Most similar substituted canonical base you will be comparing against. (default = N, any base)
can_base_train = 'N'

#Allow external validation (default false)
ext_val = False

#Balance training chunks. May be set to false for testing, otherwise set to true. 
balance_training_chunks = False

#Balance external validation chunks. May be set to false for testing, otherwise set to true. 
balance_validation_chunks = False

#Build model using Remora 
gen_model = True

#Perform model updating
model_update = False

#Number of epochs to use during training (default 50). Datatype: int
epochs = '100'

#Parameter to choose to run mod/can chunk generation 
generate_mod_chunks = False
generate_can_chunks = False

#Working Directory 
working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/CHUNK_INVESTIGATION'

#Modified chunk file, leave as '' if not generated prior 
mod_chunk_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/P2_chunks_and_models/chunk_subsets/B_0-4999/mod_chunk_merged.npz'

#Canonical chunk file, leave as '' if not generated prior
can_chunk_file = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/P2_chunks_and_models/240617_N44_p2_chunks/N_5000_pod5/chunks/can_chunk_merged.npz'

#Modified chunk validation set 
mod_val_chunk_file = ''

#Canonical chunk validation set
can_val_chunk_file = ''

#File pathway to the model checkpoint you wish to update 
model_checkpoint = ''

############################################################
#Pytorch check 

print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"CUDA version: {torch.version.cuda}")
    print(f"Number of GPUs: {torch.cuda.device_count()}")
    print(f"GPU Name: {torch.cuda.get_device_name(0)}")
else: 
    print('CUDA not compatitible with PyTorch, please verify PyTorch and CUDA drivers to make sure they are compaitible') 
    sys.exit()

#Step 0 - Initialize
print('Xemora [Status] - Initializing Xemora training...')

#Generate directories
working_dir = check_make_dir(working_dir)
chunk_dir = check_make_dir(os.path.join(working_dir, 'chunks'))
model_dir = check_make_dir(os.path.join(working_dir, 'model'))

class XemoraError(Exception):
    pass
    
def gen_training_chunks(chunk_dir, mod_chunk, can_chunk, balance_chunks): 
    """
    gen_training_chunks is a function that is merges a modified and canonical set of 
    chunks without balancing. 
    
    Parameters:
    chunk_dir: directory to output the merged training chunks 
    mod_chunk: file pathway for modified chunk (merged) 
    can_chunk: file pathway for canonical chunk (merged)
    balance_chunks: parameter to determine if the modified and canonical labels should be balanced 
    
    Return: 
    training_chunks - file pathway to generated training chunks 
    """
    
    print('Xemora [STATUS] - Merging modified and canonical chunks for training.')
    
    training_chunks = os.path.join(chunk_dir, 'training_chunks.npz')
    
    if os.path.exists(mod_chunk) and os.path.exists(can_chunk):
        if balance_chunks:
            cmd = 'remora \
              dataset merge \
              --balance \
              --input-dataset '+mod_chunk+' '+chunk_num+'_000 \
              --input-dataset '+can_chunk+' '+chunk_num+'_000 \
              --output-dataset '+training_chunks
            os.system(cmd)
        else: 
            cmd = 'remora \
              dataset merge \
              --input-dataset '+mod_chunk+' '+chunk_num+'_000 \
              --input-dataset '+can_chunk+' '+chunk_num+'_000 \
              --output-dataset '+training_chunks
            os.system(cmd)
    else:
        raise XemoraError("Modified and/or canonical chunk path is invalid, please verify file pathways")
    
    # Verify training chunks were actually generated
    if os.path.exists(training_chunks):
        print(f'Xemora [STATUS] - Training chunks successfully generated at {training_chunks}')
    else:
        print('Xemora [ERROR] - Failed to generate training chunks.')
        training_chunks = None
        sys.exit()
        
    return training_chunks 
    
def gen_validation_chunks(chunk_dir, mod_chunk, can_chunk, balance_chunks):
    """
    gen_validations_chunks is a function that is merges a modified and canonical set of 
    chunks with balancing. This will be used as the external dataset for 
    unbalancd model training
    
    Parameters:
    chunk_dir: directory to output the merged training chunks 
    mod_chunk: file pathway for modified chunk (merged) 
    can_chunk: file pathway for canonical chunk (merged)
    balance_chunks: parameter to determine if the modified and canonical labels should be balanced 
    
    Return: 
    validation_chunks - file pathway to generated validation chunks 
    """
    
    print('Xemora [STATUS] - Merging modified and canonical chunks for validation.')
    
    validation_chunks = os.path.join(chunk_dir, 'validation_chunks.npz')

    if os.path.exists(mod_chunk) and os.path.exists(can_chunk):
        if balance_chunks:
            cmd = 'remora \
              dataset merge \
              --balance \
              --input-dataset '+mod_chunk+' '+chunk_num+'_000 \
              --input-dataset '+can_chunk+' '+chunk_num+'_000 \
              --output-dataset '+validation_chunks
            os.system(cmd)
        else:
            cmd = 'remora \
              dataset merge \
              --input-dataset '+mod_chunk+' '+chunk_num+'_000 \
              --input-dataset '+can_chunk+' '+chunk_num+'_000 \
              --output-dataset '+validation_chunks
            os.system(cmd)
    else:
        raise XemoraError("Modified and/or canonical chunk path is invalid, please verify file pathways")
    # Verify training chunks were actually generated
    if os.path.exists(validation_chunks):
        print(f'Xemora [STATUS] - Validation chunks successfully generated at {validation_chunks}')
    else:
        print('Xemora [ERROR] - Failed to generate validation chunks.')
        validation_chunks = None
        sys.exit()
        
    return validation_chunks 
    
def train_model(model_dir, training_chunk_path):
    """
    train_model merges a pair of modified and canonical chunks then generates a 
    Xemora model using Remora API 
    
    Parameters: 
    model_dir: output directory for model
    training_chunk_path:
    
    Return:
    model_path: file pathway to generated model file 
    """
    
    print('Xemora [STATUS] - Training model.')
    cmd = 'remora \
      model train \
      '+training_chunk_path+' \
      --model '+ml_model_path+' \
      --device 0 \
      --output-path '+model_dir+' \
      --overwrite \
      --kmer-context-bases '+kmer_context+' \
      --chunk-context '+chunk_context+' \
      --epochs '+epochs+' \
      --val-prop '+val_proportion+' \
      --batch-size 100 '
    print(cmd)
    os.system(cmd)
    
    model_path =os.path.join(model_dir, 'model_best.pt')
    renamed_model_path = os.path.join(model_dir, mod_base_train+can_base_train+'_model.pt')
    os.rename(model_path, renamed_model_path)
     
    return renamed_model_path

def train_model_ext(model_dir, training_chunk_path, validation_chunk_path):
    """
    train_model merges a pair of modified and canonical chunks then generates a 
    Xemora model using Remora API. This allows for an external validation set to be used 
    during training for extra comparison
    
    Parameters: 
    model_dir: output directory for model
    training_chunk_path:
    
    Return:
    model_path: file pathway to generated model file 
    """
    
    print('Xemora [STATUS] - Training model.')
    cmd = 'remora \
      model train \
      '+training_chunk_path+' \
      --model '+ml_model_path+' \
      --device 0 \
      --output-path '+model_dir+' \
      --overwrite \
      --kmer-context-bases '+kmer_context+' \
      --chunk-context '+chunk_context+' \
      --epochs '+epochs+' \
      --val-prop '+val_proportion+' \
      --batch-size 100 ' + '\
      --ext-val ' + validation_chunk_path
    os.system(cmd)

    model_path =os.path.join(model_dir, 'model_best.pt')
    renamed_model_path = os.path.join(model_dir, mod_base_train+can_base_train+'_model.pt')
    os.rename(model_path, renamed_model_path)
     
    return renamed_model_path
    
def update_model(model_dir, training_chunk_path, model_checkpoint):
    """
    update_model is an optioonal method that allows you to update a provided model 
    checkpoint with new modified and canonical training chunks using Remora api. 
    """

    print('Xemora [STATUS] - Training model.')
    cmd = 'remora \
      model train \
      '+training_chunk_path+' \
      --model '+ml_model_path+' \
      --device 0 \
      --output-path '+model_dir+' \
      --overwrite \
      --kmer-context-bases '+kmer_context+' \
      --chunk-context '+chunk_context+' \
      --val-prop '+val_proportion+' \
      --batch-size 100 \
      --finetune-path '+model_checkpoint
    print(cmd)
    os.system(cmd)

    model_path =os.path.join(model_dir, 'model_best.pt')
    renamed_model_path = os.path.join(model_dir, mod_base_train+can_base_train+'_model_updated.pt')
    os.rename(model_path, renamed_model_path)
    
    return renamed_model_path

def main(): 

    if mod_chunk_file == '':
        print('Xemora [STATUS] - Modified chunk file not inputted')
        
    if can_chunk_file == '': 
        print('Xemora [STATUS] - Canonical chunk file not inputted')
        
    if gen_model and not model_update: 
        if not ext_val:
            training_chunks = gen_training_chunks(chunk_dir, mod_chunk_file, can_chunk_file, balance_training_chunks)
            output_model = train_model(model_dir, training_chunks)
            print('Xemora [Complete] - Model generated at', output_model)
        elif ext_val:
            training_chunks = gen_training_chunks(chunk_dir, mod_chunk_file, can_chunk_file, balance_training_chunks)
            validation_chunks = gen_validation_chunks(chunk_dir, mod_val_chunk_file, can_val_chunk_file, balance_validation_chunks)
            output_model = train_model_ext(model_dir, training_chunks, validation_chunks)
            print('Xemora [Complete] - Model generated at', output_model)
        else:
            print('Xemora [Error] - Please set ext_val to True or False')
            sys.exit()
    elif model_update and not gen_model: 
        '''
        this needs to be updated and ideally this is all done in 1 script instead of 2 
        '''
        training_chunks = gen_training_chunks(chunk_dir, mod_chunk_file, can_chunk_file, balance_chunks)
        output_model = update_model(model_dir, training_chunks, model_checkpoint)
        print('Xemora [Complete] - Model updated, outputted at', output_model)
    elif gen_model and model_update: 
        raise XemoraError("Cannot set both gen_model and model_update to True, please check your parameters in xr_params.py")
    else:
        print('Xemora [STATUS] - Skipping chunk merging and model training')
        sys.exist()
    
if __name__ == '__main__':
    main()
