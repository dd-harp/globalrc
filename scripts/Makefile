rc_test:
	Rscript rc_kappa_run.R --config=rc_kappa.toml --country=gmb \
	  --outvars=test --overwrite --years=2010:2011 --cores=1 --draws=3

rc_uga:
	Rscript rc_kappa_run.R --config=rc_kappa.toml --country=uga \
	  --outvars=test --overwrite --years=2010:2011 --cores=4 --draws=10

rc_world:
	Rscript rc_kappa_run.R --config=rc_kappa.toml \
	  --outvars=test --overwrite --years=2010:2011 --cores=4 --draws=3

rc_cluster_test:
	/ihme/singularity-images/rstudio/shells/execRscript.sh -i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img -s rc_kappa_run.R --config=rc_kappa.toml --country=uga --years=2010:2011  --outvars=test --overwrite --cores=24 --draws=3

rc_cluster_world_test:
	/ihme/singularity-images/rstudio/shells/execRscript.sh -i /ihme/singularity-images/rstudio/ihme_rstudio_4030.img -s rc_kappa_run.R --config=rc_kappa.toml   --outvars=test --overwrite --years=2010:2011 --cores=1 --draws=3
