"""
-----------------------------------------------------------------------------------------
post_frmiprep.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run a pre-processing analysis of pRF_gazeMod
-----------------------------------------------------------------------------------------
Prerequisite(s):
python packages: nipype, spynoza, numpy, scipy, json
imaging sofwares: fsl, freesurfer
-----------------------------------------------------------------------------------------
Example:
>>> cd /home/szinte/projects/pRF_analysis/
>>> python post_fmriprep.py "sub-002"
-----------------------------------------------------------------------------------------
"""

# IMPORTATION
# ===========
# general import
import os
import sys
import os.path as op
import glob
import json
from IPython import embed as shell
import datetime

# Nipype import
from nipype import config, logging

# Sub-workflow import
import pRF_gazeMod
from pRF_gazeMod.utils.utils import set_pycortex_config_file
from pRF_gazeMod.workflows.post_fmriprep_wf import create_post_fmriprep_workflow

# get clock
now = datetime.datetime.now()

# analysis settings
with open('analysis_settings.json') as f:
    json_s = f.read()
    analysis_info=json.loads(json_s)

os.environ['SUBJECTS_DIR'] = analysis_info['FS_subject_dir']

# change cortex database folder
pycortex_folder 	= 	os.path.join(analysis_info['FS_subject_dir'], 'cortex')
set_pycortex_config_file(project_folder 	= 	pycortex_folder)

# nipype workflow settings 
log_dir             = 	op.join(analysis_info['aeneas_project_directory'], 'log')                       # Where to store logs

config.update_config({'logging': {'log_to_file':                    True,                               # Indicates whether logging should also send the output to a file 
                                  'log_directory':                  log_dir,                            # Where to store logs
                                  'workflow_level':                 'INFO',                             # How detailed the logs regarding workflow should be 
                                  'interface_level':                'DEBUG'},                           # How detailed the logs regarding interface execution should be
                     'execution': {'stop_on_first_crash':           True}})                             # Should the workflow stop upon first node crashing or try to execute as many nodes as possible? 

config.enable_debug_mode()                                                                              # Enable debug mode
logging.update_logging(config)                                                                          # Update nipype

# subject 
sub_id              =   str(sys.argv[1])                                                                # The subject id is a commandline arguments to this script.
register_cortex     =   True                                                                            # Create the registration to pycortex

# PERFORM PRE-PROCESSING
# ======================
# Load workflow
# ------------- 
post_fmriprep_workflow = create_post_fmriprep_workflow()            									# load pre-processing workflow


# Workflow inputs
# ---------------
post_fmriprep_workflow.inputs.inputspec.fmriprep_directory              =   os.path.join(analysis_info['aeneas_project_directory'], 'fmriprep')
post_fmriprep_workflow.inputs.inputspec.output_directory                =   os.path.join(analysis_info['aeneas_project_directory'], 'pp', sub_id)
post_fmriprep_workflow.inputs.inputspec.FS_subject_dir                  =   analysis_info['FS_subject_dir'] 
post_fmriprep_workflow.inputs.inputspec.av_func                         =   analysis_info['av_across_runs_func'] 
post_fmriprep_workflow.inputs.inputspec.sub_id                          =   sub_id
post_fmriprep_workflow.inputs.inputspec.pycortex_mask_types             =   analysis_info['pycortex_mask_types']
post_fmriprep_workflow.base_dir 										= 	op.join(analysis_info['aeneas_project_directory'],'tmp', sub_id,now.strftime("%Y_%m_%d_%H_%M"))
post_fmriprep_workflow.inputs.inputspec.temp_dir_cortex             	=   op.join(analysis_info['aeneas_project_directory'],'tmp', sub_id,now.strftime("%Y_%m_%d_%H_%M"),'cortex')

# Wokflow running
# ---------------
post_fmriprep_workflow.write_graph(os.path.join(analysis_info['aeneas_project_directory'], 'pp', sub_id, 'wf.png'))
post_fmriprep_workflow.run()
