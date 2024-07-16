##############################################################################
"""
xr_convert_models.py 

By: Me :>


A convience script that outputs previously made Xemora models to be compatible 
with Dorado
"""
##############################################################################

import os
import glob
import sys

from xr_tools  import *
from xr_params import *

##############################################################################
##############################Parameters #####################################


working_dir = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/ST_Train_Q9_Dorado' #output directory for the generated model 

#xemora_model_path = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/BN_balanced_training_copy/BN_model.pt' #should be a .pt file

#xemora_model_checkpoint = '/home/marchandlab/github/jay/xemora_randomer_troubleshooting/xemora_models/BN_balanced_training_copy/model_best.checkpoint' #should be a .checkpoint file
xemora_model_checkpoint = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/BSn/240110_BSn_xr_train/ST_Train_Q9/model/model_best.checkpoint'
##############################################################################
#Generate directories
print('Xemora [Status] - Initializing Xemora basecalling.')

working_dir = check_make_dir(working_dir)

##############################################################################

print('Xemora [STATUS] - Converting existing Xemora model to a dorado model')

#cmd = 'remora model export --model-path '+xemora_model_path+' '+xemora_model_checkpoint+' '+working_dir
cmd = 'remora model export '+xemora_model_checkpoint+' '+working_dir
        
print(cmd)
os.system(cmd)
