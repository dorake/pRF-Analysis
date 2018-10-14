# Analysis software for 7T multiband pRF_gazeMod experiment

This is the analysis repository of the pRF_gazeMod exeperiemnt, in which pRF are measured for different gaze posiiton.
Raw data are converted in BIDS format using bids_creator.py and will be uploaded to openfMRI.org. 

The software to run the experiment itself is stored in [https://github.com/mszinte/pRF_gazeMod_testCode]


## Analysis specifics

- Data are first pre-processed using fmriprep on LISA server by runnnig fmriprep_lisa.py.
- freesurfer segmentation is corrected manually (see segmentation.sh)
- surface must be cut for flateting (see segmentation.sh)
- pycortex database is created out of the freesurfer segmentation
- data are then processed using a custom preprocessing workflow: post_frmiprep.py
- pRF parameters are extracted using pRF_gazeMod/cartesius_lisa/prf_adapt_jobscript.py on Cartesius or Lisa.
- pRF parameters are analysed and flatmap printed with pp_roi.py
- ROI are drawn using Inkscape
- basic pRF analysis plots (pRFmap, pRFecc, pRFcov, pRFtime plots,) are produced with main_analysis.py
- eye data and behavior plots are produced with ?
- ...

## To do

- add error bars in bar plots in ROI figure
- bar plots into ROI basic
- ROI main: x_gain, x_diff, amplitude_gain

- check correct fit with og model
- do change in main_analysis.py
- markdown of cuts for flatten maps
- markdown document for ROI drawing
- general bash code runnning each piece together

