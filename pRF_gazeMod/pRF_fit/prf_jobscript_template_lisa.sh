#!/bin/bash
#PBS -lwalltime=---job_dur---
#PBS -lnodes=1
cd /home/pRF_analysis/
python -m ---fit_file--- ---slice_no--- ---subject--- ---train_input--- ---test_input--- ---path--- ---gaze_condition---
