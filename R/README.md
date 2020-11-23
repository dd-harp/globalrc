# Generate Scaled Attack Rate


## Files in this directory

* `AgeOperators.R` - makes a matrix operator to age a population.
* `AgeWeights.R` - Relative size of population by age group.
* `MetricConversions.R` - Probability-related conversions of values.
* `adam.R` - The model, from Ramp-model-library.
* `combine_pr_files.R` - a script to take the input PR data and put it in one CSV for reading.
* `gen_scaled_ar.R` - Script that scales AR. This runs as a cluster job.
* `launch_gen_scaled_ar.$` - Interactive notes on how to run the gen_scaled.
* `print_launch.R` - Prints the qsub that launches gen_scaled_ar.R.
* `run_ramp.md` - How to login to the cluster and run this stuff.
* `scale_ar.R` - functions that do the scaling by case rates. gen_scaled_ar.R calls this.
* `test-adam.R` - tests for the adam model
* `test-adam_series.R` - tests for one function off the adam model.
* `test-scale_ar.R` - tests for the scale_ar to rescale by case rates.

## The rc_kappa code.

Look at `rc_kappa.md` for more.
