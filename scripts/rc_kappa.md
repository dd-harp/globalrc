# Running rc_kappa

## Introduction

This is how to run the `rc_kappa` code. There are two ways this code runs,
in parallel on a single computer, and distributed on the cluster.

## What's in the package

### Shell Scripts

### R main scripts and shell scripts

* `rc_kappa.sh` - A shell script to qsub the single-machine version of
  the code. This was failing on the cluster because R's fork-cluster
  would get a SIGPIPE singal.
* `rc_worker.sh` - A shell script to qsub for distributed runs.
* `rc_kappa_run.R` - A main() to run the single-machine version of the code.
* `rc_run_plan.R` - A main() to run the first distributed step.
* `rc_run_worker.R` - What runs in tasks for the main distributed step.
* `rc_run_assemble.R` - A main() to run the distributed assembly step.

### Code

* `rc_kappa.R` - The main file with most of the code.
* `serialize.R` - Saves and loads HDF chunks of data, for each tile.
* `plan_io.R` - Saves and loads the plan to decompose the data for parallel
  computation.
* `hilbert.R` - Orders tiles so that nearby tiles are nearby in a linear list.
* `test-*.R` - Testing files for the code.

### Science

* `rc_kappa.Rmd` - A vignette to look at the per-pixel function that decides
  the numbers. All the science is here.
* `rc_kappa.toml` - Parameters for the run.
* `rc_sensitivity.jl` - Calculates Jacobian of the pixel work and
  does Bayesian inference on it.

## Running distributed on the cluster

There is some introductory information for using the cluster in
`run_ramp.md`. For `rc_kappa`, in particular, there are three steps.

1. Create a plan.
2. Run workers to compute that plan.
3. Assemble the results from the workers.

All three steps take the same command line, but they have different scripts.

### Create a plan.

The first script is `rc_run_plan.R`.

* `--config=rc_kappa.toml` This one uses parameters from the
  `rc_kappa.toml` configuration file.
* `--years=2000:2019` You have to give it a years argument,
  and that argument must be at least two years, for instance, 2017:2018.
* `--draws=100` The draws can be 1, in which case it uses mean values.
* `--tasks=100` Not used for this step, but I like to keep the arguments
  consistent. This doesn't match draws in any way. It's the number of
  SGE jobs we will run.
* `--outvars=201121_split100` This is the output directory. It will
  be: `/ihme/malaria_modeling/projects/globalrc/outputs/basicr/201121_split100`.

```
/ihme/singularity-images/rstudio/shells/execRscript.sh -i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img -s rc_run_plan.R --config=rc_kappa.toml --outvars=201121_split100 --years=2000:2019 --draws=100 --tasks=100
```

This step is quick to run. I do it interactively. It creates two files in
the `outvars` directory:

* `plan.json` - This is the extent of the map to process and information
  about the blocksize and coordinates in lat-long.
* `work.csv` - This is a list of tiles which have a nonzero PfPR value,
  where a tile is a 32x32 set of pixels in each image. The blocksize
  is set in the `rc_kappa.toml` file.

### Run workers on that plan

The worker runs under SGE. The shell script to qsub, with `qsub rc_worker.sh`,
is set up to handle both the initial run and rerunning any task that failed
to complete. In it, you'll see these choices.

* `-cwd` - Start the process in the same directory as the one in which
  I invoke the `qsub rc_worker.sh`.
* `fthread=2` - Each R process takes one thread, but it can be helpful
  to allow another thread because I/O can take its own thread sometimes.
* `-o /share/temp/sgeoutput...` - This says that each task gets its
  own output file, containing both output and error, because of the
  `-y` argument.
* `-t 1-100` - This tells the scheduler to run 100 copies of this script
  as what the scheduler calls tasks. Each task has its own `SGE_TASK_ID`
  set to its index in the 100, starting from 1, going to 100.
  Each time the R code runs, it checks the `SGE_TASK_ID` environment
  variable in order to figure out its task.
* `-q all.q` - This is the general queue. Works fine.
* `-l h_rt=3:00:00` - Allow three hours to run. One hour was too litte,
  and sometimes I/O can be terribly much slower than normal, so leave room,
  but not 72 hours, because they it won't get scheduled soon.
* `-l m_mem_free=4G` - Request 4GB of RAM. This seems to need about 2GB.

This worker took over an hour for 9 tiles per worker for 100 draws.
That means 1 tile for 1000 draws could be 2 hours.

If I haven't run this script on a given number of draws, I run it
first for tasks `-t 100-100`, and record its job id,
in order to see how long it takes and how much memory it needs.
These are in `qacct -j <job_id>` as RSS (in kb) and wallclock time.
Another way, easier, is to run the comman interactively, prefixed
by `/usr/bin/time --verbose` and adding `--task=1` at the end.
This will both the real time and
the maximum resident set size, in kilobytes. Divide by 1024^2 to
get GB for the qsub command. It should be near 1-4 GB.

```
#!/bin/sh
#$ -P proj_mmc
#$ -cwd
#$ -j y
#$ -o /share/temp/sgeoutput/adolgert/rc_kappa_$JOB_ID_$TASK_ID.txt
#$ -N rc_split
#$ -l fthread=2
#$ -t 1-100
#$ -q all.q
#$ -l h_rt=3:00:00
#$ -l m_mem_free=4G
#$ -S /bin/bash

# If there were tasks that didn't run, paste their numbers here.
# This way you can submit one job that works through the missing ones.
# MISSING=40,44,47,49,50,52
MISSING=
if [[ -n "${MISSING}" ]]
then
  export SGE_TASK_ID=`echo $MISSING | cut -d"," -f"${SGE_TASK_ID}"`
fi

VERSION=201121_split100
/ihme/singularity-images/rstudio/shells/execRscript.sh \
  -i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img \
  -s rc_run_worker.R --config=rc_kappa.toml --outvars=${VERSION} \
  --years=2000:2019 --tasks=100 --draws=100
```

The `MISSING` lines are there in case any tasks fail. This is a common
problem on this cluster. You'll see that the next step, to assemble
the data, will report missing tasks to the command line. You can
copy those values here. Then, if there are 6 missing, change the tasks
to go from `-t 1-6`, and, instead of executing the 1st-6th tasks,
it will do the missing ones.

These runs generate one HDF5 file for each worker, in the same
outvars directory as above.

### Assemble output from workers

Finally, we read all of the HDF5 files and create GeoTIFFs, one for
each variable, for each median and confidence interval, for each year.

The `rc_run_assemble.R` script takes the same arguments as those above.
If this script doesn't find all of its inputs, it will list those
that are missing and print the missing tasks in a format suitable
to place in the `MISSING` variable in the qsub script above.
If you know some jobs failed, running this is any easy way to find the
failures.

This script takes about 3GB to read data from 100 files.
It works by variable, so the memory usage is low. It would be possible
to run this in parallel over the variables. Takes about 20 mins now.

```
/ihme/singularity-images/rstudio/shells/execRscript.sh -i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img -s  rc_run_assemble.R --config=rc_kappa.toml --outvars=201121_split100 --years=2000:2019 --draws=100 --tasks=100
```

The workers can run in parallel with GNU parallel. You use the task
ID to say which variable this worker will save. There are currently
18 variables, but may change.

```
parallel Rscript scripts/rc_run_assemble.R --config=scripts/sam_mean.toml   --outvars=201124_africa_mean --years=2017:2017 --draws=100 --tasks=217 --overwrite --task={} ::: {1..18}
```
