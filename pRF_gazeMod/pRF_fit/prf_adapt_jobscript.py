"""
-----------------------------------------------------------------------------------------
prf_adapt_jobscript
-----------------------------------------------------------------------------------------
Goal of the script:
create jobscript to run in a cluster (LISA or CARTESIUS)
note that cd must be where this code is cd /pRF_analysis/pRF_fit/
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name (e.g. 'sub-001')
sys.argv[2]: gaze condition (e.g. 'gazeRight')
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
for s in sub-001 sub-002 sub-003 sub-004
do
    cd /home/projects/pRF_gazeMod/pRF_analysis/pRF_fit/
    python prf_adapt_jobscript.py $s gazeRight
    python prf_adapt_jobscript.py $s gazeLeft
done
-----------------------------------------------------------------------------------------
"""

# General imports
# ---------------
import numpy as np
import re
import os
import glob
import json
import sys
import pdb
import datetime

# MRI import
# ----------
import nibabel as nb

# Initial settings
# ----------------
# load the analysis parameters from json file
with open(os.path.abspath('../../analysis_settings.json')) as f:                                            # get main analysis settings
    json_s                      =   f.read()
    analysis_info               =   json.loads(json_s)

subjects                        =   [sys.argv[1]]                                                           # subject name (e.g. sub-001)
gaze_condition                  =   sys.argv[2]                                                             # gaze condition name (e.g gazeRight)

if analysis_info['server_name'] == 'cartesius':
    jobscript_template_file     =   os.path.join(os.getcwd(), 'prf_jobscript_template_cartesius.sh')        # jobscript to replace with parameters
    basedir                     =   analysis_info['cartesius_project_directory']                            # output directory
    fitperhour                  =   1250.0                                                                  # number of fit per hours that cartesius can do
elif analysis_info['server_name'] == 'lisa':
    jobscript_template_file     =   os.path.join(os.getcwd(), 'prf_jobscript_template_lisa.sh')             # jobscript to replace with parameters
    basedir                     =   analysis_info['lisa_project_directory']                                 # output directory
    fitperhour                  =   1000.0                                                                  # number of fit per hours that lisa can do
elif analysis_info['server_name'] == 'aeneas':
    jobscript_template_file     =   os.path.join(os.getcwd(), 'prf_jobscript_template_aeneas.sh')           # jobscript to replace with parameters
    basedir                     =   analysis_info['aeneas_project_directory']                               # output directory
    fitperhour                  =   500.0                                                                   # number of fit per hours that lisa can do

fit_script                      =   'pRF_gazeMod.pRF_fit.prf_one_fold'                                      # script to run

# Create job folders
# ------------------
try:
    os.makedirs(                    os.path.join(basedir, 'pp', 'jobs'))
except:
    pass
os.chdir(                           os.path.join(basedir, 'pp', 'jobs'))

# Main loop
# ---------
for subject in subjects:
    # Determine slice number from cortical mask file
    mask_data                   =   nb.load(os.path.join(basedir, 'pp', subject, 'masks',                   # get the cortical mask
                                            'cortex_cortical.nii.gz')).get_data()
    slices                      =   np.arange(mask_data.shape[2])[mask_data.mean(axis=(0,1))>0]             # get only slice with voxels

    # Define training and test inputs
    train_inputs                =   sorted(glob.glob(os.path.join(basedir, 'pp', subject,                   # list of trainning input: leave one out averages 
                                            'av_' + gaze_condition, 'loo', '*.nii.gz')))

    test_inputs                 =   sorted(glob.glob(os.path.join(basedir, 'pp', subject, 'tf', 'Z',        # single run test input
                                            '*%s*.nii.gz'%gaze_condition)))

    # Loop across trainning and test inputs
    for (train_input, test_input) in zip(train_inputs, test_inputs):
        for slice_no in slices:
            jobscript           =   open(jobscript_template_file)                                           # open template script to edit
            working_string      =   jobscript.read()
            jobscript.close()
            num_vox             =   mask_data[:, :, slice_no].sum()                                         # get the number of voxels in slice
            job_dur             =   str(datetime.timedelta(hours=np.ceil(num_vox/fitperhour)))              # compute estimated job duration

            base_file_name      =   os.path.split(train_input)[-1][:-7]                                     # define base_file_name
            opfn                =   os.path.join(basedir,'pp',subject, 'cv_'+gaze_condition, 'prf',         # define output file name
                                                        base_file_name + '_est_%s.nii.gz' % str(slice_no).zfill(2))
            if os.path.isfile(opfn):
                if os.path.getsize(opfn) != 0:
                    print('output file %s already exists and is non-empty. aborting slice %i.'%(opfn,slice_no)) # quit code if data exist
                    continue

            RE_dict = {             '---job_dur---':            job_dur,                                    # duration of the job
                                    '---fit_file---':           fit_script,                                 # script to use
                                    '---subject---':            subject,                                    # subject name
                                    '---slice_no---':           str(slice_no),                              # slice number
                                    '---train_input---':        train_input,                                # trainning data path
                                    '---test_input---':         test_input,                                 # test data path
                                    '---path---':               basedir,                                    # output data path
                                    '---gaze_condition---':     gaze_condition}                             # condition name

            for e in RE_dict.keys():
                working_string  =   working_string.replace(e, RE_dict[e])                                   # replace sting in jobscript

            js_name             =   os.path.join(basedir, 'pp', 'jobs', subject + '_' + str(slice_no) +     # jobscript bash filename (8:-7 to get rid of subject name and .nii.gz)
                                            '_' + os.path.split(train_input)[-1][8:-7].replace('/', '-') + 
                                            '.sh')
            
            of                  =   open(js_name, 'w')                                                      # create bash file
            of.write(working_string)                                                                        # write command  
            of.close()                                                                                      # close bash file
            
            print(                  'submitting ' + js_name + ' to queue')                                  # print message
            
            if analysis_info['server_name'] == 'lisa':                                                      
                os.system('qsub ' + js_name)                                                                # run bash file in lisa
            elif analysis_info['server_name'] == 'cartesius':
                os.system('sbatch ' + js_name)                                                              # run bash file in cartesius
            elif analysis_info['server_name'] == 'aeneas':
                os.system('sh js_name')

