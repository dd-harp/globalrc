#!/bin/sh
#$ -P proj_mmc
#$ -cwd
#$ -j y
#$ -o /share/temp/sgeoutput/adolgert/rc_kappa_$JOB_ID_$TASK_ID.txt
#$ -N rc_split
#$ -l fthread=2
#$ -t 1-2
#$ -q all.q
#$ -l h_rt=3:00:00
#$ -l m_mem_free=4G
#$ -S /bin/bash

# If there were tasks that didn't run, paste their numbers here.
# This way you can submit one job that works through the missing ones.
# MISSING=40,44,47,49,50,52
MISSING=134,212
if [[ -n "${MISSING}" ]]
then
  export SGE_TASK_ID=`echo $MISSING | cut -d"," -f"${SGE_TASK_ID}"`
fi

VERSION=201122_1000
/ihme/singularity-images/rstudio/shells/execRscript.sh \
  -i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img \
  -s  rc_run_worker.R --config=rc_kappa.toml --outvars=${VERSION} \
  --years=2000:2019 --draws=1000 --tasks=1010
