rc_test:
	Rscript rc_kappa_run.R --config=rc_kappa.toml --country=gmb \
	  --outvars=test --overwrite --years=2010:2011 --cores=1 --draws=3

rc_single:
	Rscript rc_kappa_run.R --config=rc_kappa.toml --country=gmb \
	  --outvars=test --overwrite --years=2011:2011 --cores=1 --draws=3

rc_uga:
	Rscript rc_kappa_run.R --config=rc_kappa.toml --country=uga \
	  --outvars=test --overwrite --years=2010:2011 --cores=4 --draws=10

rc_uga_single:
	Rscript rc_kappa_run.R --config=rc_kappa.toml --country=uga \
	  --outvars=test --overwrite --years=2010:2011 --cores=4 --draws=10

rc_world:
	Rscript rc_kappa_run.R --config=rc_kappa.toml \
	  --outvars=test --overwrite --years=2010:2011 --cores=4 --draws=3

rc_cluster_test:
	/ihme/singularity-images/rstudio/shells/execRscript.sh \
	-i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img \
	-s rc_kappa_run.R --config=rc_kappa.toml --country=uga \
	--years=2010:2011  --outvars=test --overwrite --cores=24 --draws=3

rc_cluster_world_test:
	/ihme/singularity-images/rstudio/shells/execRscript.sh \
	-i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img \
	-s rc_kappa_run.R --config=rc_kappa.toml   --outvars=test \
	--overwrite --years=2010:2011 --cores=1 --draws=3

sam_mean:
	Rscript rc_kappa_run.R --config=sam_mean.toml \
	  --outvars=201123_sam_mean --years=2017:2017 --cores=4 --draws=3

world_mean:
	Rscript rc_kappa_run.R --config=world_mean.toml \
	  --outvars=201123_global_median --years=2019:2019 --cores=4 --draws=3

world_mean_plan:
	/ihme/singularity-images/rstudio/shells/execRscript.sh \
	-i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img \
	-s rc_run_plan.R --config=world_mean.toml   --outvars=201123_global_median --years=2019:2019 --draws=100

world_mean_worker:
	/ihme/singularity-images/rstudio/shells/execRscript.sh \
	-i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img \
	-s rc_run_worker.R --config=world_mean.toml   --outvars=201123_global_median \
	--years=2019:2019 --draws=3 --tasks=9000

render:
	rmarkdown::render("notebooks/check_pixel_work.Rmd")
	