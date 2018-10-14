# Use rsync to copy data from/to server
# rsync -a --no-g --no-p -vzhe  -o Compression=no -x --progress server:/home/ /home

import re
import os
import glob
import json
import sys

subjects = ['001','002','003','004']

batch_string = """# shell for the job:
#PBS -S /bin/bash
#PBS -lwalltime=50:00:00 -lnodes=1:mem64gb
# job requires at most 50 hours, 0 minutes
#     and 0 seconds wallclock time

# call the programs
echo "Job $PBS_JOBID started at `date`" | mail $USER -s "Job $PBS_JOBID"

PYTHONPATH="" singularity run /home/poldracklab_fmriprep_1.0.0-2017-12-07-ad31fad01869.img /home/sourcedata /home/derivatives participant --participant-label SJ_NR -w $TMPDIR/pRF_gazeMod/work --output-space T1w fsaverage5 --nthreads 15 --use-syn-sdc --low-mem --fs-license-file /home/software/freesurfer/license.txt

wait          # wait until programs are finished

echo "Job $PBS_JOBID finished at `date`" | mail $USER -s "Job $PBS_JOBID"
"""

basedir = '/home/batch/'

os.chdir(basedir)


for subject in subjects:

    working_string = batch_string.replace('SJ_NR', subject)

    js_name = os.path.join(basedir, subject + '.sh')
    of = open(js_name, 'w')
    of.write(working_string)
    of.close()

    print('submitting ' + js_name + ' to queue')
    os.system('qsub ' + js_name)

