# Population Receptive Field Analysis

Scripts for fMRI (pre-)processing and analysis of a population receptive field (pRF) experiment that I worked on as part of my master thesis and internship at the Spinoza Centre for Neuroimaging, Amsterdam, Netherlands. Scripts assume a config .json file and the data must be in standardized BIDS format (http://bids.neuroimaging.io). Raw data can be converted to BIDS using bidsCreator.py


## Analysis specifics

- Data are first pre-processed using fmriprep (https://github.com/poldracklab/fmriprep) on LISA server (https://userinfo.surfsara.nl/systems/lisa) by runnnig fmriprep_lisa.py.
- freesurfer segmentation is corrected manually
- surface must be cut for flattening
- pycortex database is created out of the freesurfer segmentation
- data are then processed using a custom preprocessing workflow: post_frmiprep.py
- pRF parameters are extracted using pRF_gazeMod/cartesius_lisa/prf_adapt_jobscript.py on Lisa.
- pRF parameters are analysed and flatmap printed with pp_roi.py
- ROI are drawn using Inkscape
- basic pRF analysis plots (pRFmap, pRFecc, pRFcov, pRFtime plots,) are produced with main_analysis.py
- eye data plots are produced with eye_analysis.ipynb

