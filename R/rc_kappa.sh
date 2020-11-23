#!/bin/sh
#$ -P proj_mmc
#$ -cwd
#$ -j y
#$ -o /share/temp/sgeoutput/adolgert/rc_kappa_$JOB_ID.txt
#$ -N rc_short
#$ -l fthread=24
#$ -q all.q
#$ -l h_rt=2:00:00
#$ -l m_mem_free=100G
#$ -S /bin/bash

VERSION=`date +%y%m%d_%H%M%S`_${JOB_NAME}
CORES=${SGE_HGR_fthread}
/ihme/singularity-images/rstudio/shells/execRscript.sh -i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img -s rc_kappa_run.R --config=rc_kappa.toml --outvars=${VERSION} --years=2000:2019 --cores=${CORES} --draws=100
