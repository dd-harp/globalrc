#!/bin/sh
#$ -P proj_mmc
#$ -cwd
#$ -j y
#$ -o /share/temp/sgeoutput/adolgert/rc_kappa_$JOB_ID_$TASK_ID.txt
#$ -N rc_global
#$ -l fthread=2
#$ -t 1-${tasks}
#$ -q all.q
#$ -l h_rt=3:00:00
#$ -l m_mem_free=4G
#$ -S /bin/bash
#
# If there were tasks that didn't run, paste their numbers here.
# This way you can submit one job that works through the missing ones.
# MISSING=40,44,47,49,50,52
MISSING=
if [[ -n "${MISSING}" ]]
then
  export SGE_TASK_ID=`echo $MISSING | cut -d"," -f"${SGE_TASK_ID}"`
fi
/ihme/singularity-images/rstudio/shells/execRscript.sh \
	-i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img \
	-s rc_run_worker.R --config=${config} \
	--outvars=${outvars} \
	--years=${years} 	${country} \
	--draws=${draws} --tasks=${tasks}
