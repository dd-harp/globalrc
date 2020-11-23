# Run RAMP Analytics on the Cluster

In order to run our code on the cluster, you have to get on the 
cluster, configure R, get the code, and work with the scheduler.
Here are the steps.

## Get on the cluster

1. Connect to the cluster using Big-EDGE IP VPN.
2. SSH to the cluster.
```
ssh cluster-submit2.ihme.uw.edu
```

3. Login to a node. Type this. Wait. You will then be on a node where you can work.
```
qlogin -l h_rt=12:00:00,fthread=2,m_mem_free=4G -q i.q -P proj_mmc -now no
```

If you have one qlogin session and want a second one, look up the hostname
bu running the command `hostname` in the qlogin,
do another ssh to the `cluster-submit1` or `cluster-submit2`, and ssh with
```
ssh int-uge-archive-p004.cluster.ihme.uw.edu
```
I'm telly you this because it's hard to figure out the full hostname has
five parts.

## Set up R
1. Put the code where Austin puts the code.
```
mkdir -p dev
cd dev
git clone git@github.com:dd-harp/analytics-pipeline.git
git clone git@github.com:dd-harp/pr2ar.git
```

Check that R is cofigured. You configure R on the cluster by installing
an Rprofile, like this one.
```
$ cat ~.Rprofile
.libPaths(c("~/share/R3.6.2", .libPaths()))
options(repos = structure(c(CRAN = "https://repo.miserver.it.umich.edu/cran")))
```
That tells R where to put its libraries when it runs.

Copy a script to run R. IHME provides a command to run R,
but it doesn't read the Rprofile. We hacked it to load the .Rprofile so that
we can a) install our own external packages and b) use our code as packages.
The three parts are

- `~adolgert/bin/rsing` - command-line R
- `~adolgert/bin/execRscript.sh` - This runs a script as `execRscript.sh <scriptname>`.
- `~adolgert/bin/rsingbash` - A Bash command line within the R singularity
  image. This helps you figuure out what your script sees when it runs.
  For instance, `git` isn't installed on the singularity image, so you 
  can't use `git`.

So copy those if you need them.

## Install packages

```
devtools::install("~/dev/pr2ar")
```
If you're re-installing pr2ar a million times because you're editing it,
try an alternative installation call.
```
devtools::install("~/dev/pr2ar", dependencies = FALSE, upgrade = FALSE)
```

If you see compilation messages like this, everything looks fine, but you might be in trouble.
```
g++  -std=gnu++11 -I"/usr/local/lib/R/include" -DNDEBUG -I. -I/usr/local/include   -UDEBUG -DNDEBUG -DU_HAVE_ELF_H  -w -fno-gnu-unique -fno-optimize-sibling-calls -DMKL_ILP64 -m64 -I/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/include -I/usr/include                             -msse2 -mfpmath=sse -g -O3  -fPIC -fPIC   -w -fno-gnu-unique -fno-optimize-sibling-calls -DMKL_ILP64 -m64 -I/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/include -I/usr/include                             -msse2 -mfpmath=sse -g -O3   -c stri_trans_normalization.cpp -o stri_trans_normalization.o
```
The problem is that the sse and sse2 optimizations may work on one computer but not another. If you later see an error called "illegal instruction," then you have to go back and install every package with optimization turned off. This means that R's install won't work by default. You need to customize installation for failing packages. And maybe tell Infra that you'd like a cluster that uses the same processors for all of the nodes.

Test that it loads with `library(pr2ar)`.


## Run RStudio

Go the Hub and search for "R/RStudio Cluster Cheat Sheet".
Or try this link:
https://hub.ihme.washington.edu/pages/viewpage.action?pageId=60485097

Then run
```
singRstudio -P proj_mmc
```
You'll see something like this, where the last line tells you where to
point your web browser.
```
$ singRstudio -P proj_mmc
--------------------------------------------------------------------------------
RStudio Qsub Script Info
----------------------------------------
Script Name             : 'rstudio_qsub_script'
Version                 : 2.0.1

Number of threads set to: 2 (default)
Maximum memory set to   : '1G' (default)
Queue set to            : 'all.q' (default)
Runtime set to          : '03:00:00:00' (default)
Cluster set to          : 'prod' (default)
Project set to          : 'proj_mmc'
Port set to             : '7285' (based on user name)
Image type              : LBD Team (default)
Image version set to    : /share/singularity-images/lbd/rstudio/shells/lbd_rstudio_shells/rstudio_shell_9999_current.sh
Requesting archive nodes: TRUE (default)

Submitting RStudio job with the following command:
'qsub -N rstudio_ide -j y -o /share/temp/sgeoutput/adolgert/output -now no -P proj_mmc -l fthread=2 -l m_mem_free=1G -l h_rt=03:00:00:00 -q all.q  -l archive=TRUE  /share/singularity-images/lbd/rstudio/shells/lbd_rstudio_shells/rstudio_shell_9999_current.sh 7285'
Your job 40923023 ("rstudio_ide") has been submitted
http://gen-uge-archive-p207.cluster.ihme.washington.edu:7285
```
That last line is a URL. Point your web browser there, and you'll see R.

Is it pointing to your installed libraries? Check for your library with
```
.libPaths()
```
I had a problem with the versions not matching, so I ran this to make them match.
```
.libPaths(c("/ihme/homes/adolgert/share/R3.6.1", "/usr/local/R-3.6.1/library"))
```

I finally found the post that leads to how you run Rstudio on the cluster
using the Scicomp version:

	Hi, RStudio singularity users, we just released 3.6.0.2 with bug fixes and new packages. Here is the release note: https://hub.ihme.washington.edu/pages/viewpage.action?spaceKey=DataScience&title=GBD+Rstudio+Singularity+Image+Release+Note. Since it's still in beta, the default image is still 3.6.0.1. If you want to use this image, please specify -i /ihme/singularity-images/rstudio/ihme_rstudio_3602-beta.img on your command. Example:

	sh /ihme/singularity-images/rstudio/shells/jpy_rstudio_qsub_script.sh -i /ihme/singularity-images/rstudio/ihme_rstudio_3602-beta.img -t rstudio -P proj_burdenator -G w

Maybe that helps. The previous instructions are for using the geospatial version
because they have instructions.


## Fixing Permission Problems

Sometimes the permissions aren't right. These are commands that fix them.
Note that the ihme-malaria group id is 700187

```
cd /ihme/malaria_modeling
find . -type d -user adolgert -group "Domain Users" -exec chown adolgert:ihme-malaria {} \;
find . -type f -user adolgert -group "Domain Users" -exec chown adolgert:ihme-malaria {} \;
find . -user adolgert -exec chmod ug+rwX {} \;
```


## Running jobs
Modify the `out_dir` in the gen_scaled_ar.R script so that it writes to a new location.

Get the qsub string by running a script to print it.
```
execRscript.sh print_launch.R
```
Then copy the qsub string and paste it into the command line. It will be something like this.
```
qsub -l m_mem_free=2.0G -l fthread=1 -l h_rt=24:00:00 -l archive=True -q all.q -cwd -P proj_mmc -e /share/temp/sgeoutput/adolgert/errors -o /share/temp/sgeoutput/adolgert/output -N gen_scaled_ar -t 1:13500 /homes/adolgert/bin/execRscript.sh /homes/adolgert//dev/analytics-pipeline/gen_scaled_ar/gen_scaled_ar.R
```
It isn't running from inside R because `SGE_ROOT` is defined as `/opt/sge` which isn't in the
singularity image. Whatever. Run the command from the command line.
Look at the job ID that's returned. Use `qstat` to see the jobs running. When they are done,
use `qacct -j <jobid> | grep exit_status | uniq` to see if they are all exit status 0.

The files should be in `out_dir` from gen_scaled_ar.R, which is currently set to
`/ihme/malaria_modeling/projects/uganda2020/outputs/practice_draws/`.
See that they were made and have the right owner and group and privileges,
like this, so ihme-malaria owns them and they are rw for owner and group.
```
-rw-rw-r-- 1 adolgert ihme-malaria 6619 May  2 18:13 1_1.csv
```

Combine the output results by running another script. Edit this script first so that
it knows where to put its outputs without overwriting.
```
execRscript.sh combine_pr_files.R
```
Do you want to know how long it took to run jobs?
After they are finished, run
```
qacct -j gen_scaled_ar | grep ^wallclock | cut -d' ' -f5
```
This gives seconds for each run. It is about 2 hours currently, per draw.


## Data and Plots

The combined files are in a directory on the cluster, determined by the
`combine_pr_files.R` script. You can leave them there
or zip and transfer them, but go into the `plot` subdirectory to find
`make_maps.Rmd` in order to make maps from those plots. You'll point
a variable there at the location of the files.


## RStudio

If you want to run under Rstudio on the cluster, then create a ~/.Rprofile
script and `export R_PROFILE=$HOME/.Rprofile` in your .bashrc.
That script should point to a library path that matches the current
version of R in R singularity. That's `~/share/R3.6.3` right now.

Run Rstudio with this qsub for 2 threads, 10GB of memory 12 hours, interactive queue.
Yes, they changed the specifications from the way qsub handles them.
Read the instructions and then try a bunch of times until it works.
```
/ihme/singularity-images/rstudio/shells/jpy_rstudio_qsub_script.sh -t rstudio -i /ihme/singularity-images/rstudio/ihme_rstudio_3630.img -f 2 -o 2 -m 10G -h 12 -P proj_mmc -q i -l /share/temp/sgeoutput/$USER/errors
```
When you run qsub, it will tell you the port number it will try to give you,
although that can fail, of course. But then you qstat, which might look
like this:
```
$ qstat
job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID 
------------------------------------------------------------------------------------------------------------------------------------------------
  41622257 0.50026 QLOGIN     adolgert     r     05/10/2020 16:46:07 i.q@int-uge-archive-p005.clust                                    1        
  41624723 0.50023 rst_ide_20 adolgert     r     05/10/2020 17:12:22 i.q@int-uge-archive-p004.clust                                    1        

```
That tells me to connect my web browser to int-uge-archive-p004.cluster.ihme.washington.edu:6232,
because that's the port number associated with my name. It may be someone else's too, but whatever.
