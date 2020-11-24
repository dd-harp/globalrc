
library(futile.logger)


#' Asks R for a stack trace that shows which function had the problem.
#' 
#' https://renkun.me/2020/03/31/a-simple-way-to-show-stack-trace-on-error-in-r/
#' @export
improved_errors <- function() {
    options(error = function() {
        sink(stderr())
        on.exit(sink(NULL))
        traceback(3, max.lines = 1L)
        if (!interactive()) {
            q(status = 1)
        }
    })
}


#' Given a list of matrices, construct a single three-dimensional array.
#' @param list_of_matrices Each item in the list is a two-dimensional array of the same size.
#' @return A single three-dimensional array where the first dimension is the
#'     number of items in the list and the next two are the dimensions of each array.
collect_and_permute <- function(list_of_matrices) {
    sentinel <- dim(list_of_matrices[[1]])
    together <- array(do.call(c, list_of_matrices), dim = c(sentinel, length(list_of_matrices)))
    aperm(together, c(3, 1, 2))
}


#' Given the raw input data, prepare it for computation.
#' @param data A list where each item is input data.
#' @param work The matrix of which tiles to do.
#' @param process_extent Bounding box of loaded data within the computational domain.
#' @param blocksize The size of each tile.
#' @return A list with transformed data items. This list is broken
#'     up into the work tiles, so that we don't transfer the whole
#'     image to each processor when we do calculations.
#'
#' This is data movement to set up for the next algorithm.
prepare_timeseries <- function(data, work, process_extent, blocksize) {
    stopifnot(is_tile(blocksize))
    pfpr <- collect_and_permute(data$pfpr)
    am <- collect_and_permute(data$am)
    stopifnot(length(dim(pfpr)) == 3)
    stopifnot(length(dim(am)) == 3)

    work_cnt <- dim(work)[2]
    pieces <- lapply(1:work_cnt, FUN = function(col) {
        tile <- work[, col]
        tile_domain_extent <- tile_extent(tile, blocksize)
        # relative extent to the subset of data that is loaded.
        rel <- find_relative_extent(process_extent, tile_domain_extent)
        # flog.debug(paste("tile", paste(tile, collapse = ","),
        #     "relative extent", paste(rel, collapse = ","), "\n"
        # ))
        pfpr_chunk <- pfpr[, rel["rmin"]:rel["rmax"], rel["cmin"]:rel["cmax"], drop = FALSE]
        am_chunk <- am[, rel["rmin"]:rel["rmax"], rel["cmin"]:rel["cmax"], drop = FALSE]
        list(tile = tile, pfpr = pfpr_chunk, am = am_chunk)
    })

    list(
        parameters = data$parameters,
        pr_to_ar_dt = data$pr_to_ar_dt,
        years = data$years,
        chunks = pieces
        )
}


#' If you want to run from the R command line, you can call this.
#' @param args Command-line arguments, after they've been checked.
#' 
#' Use the same arguments as the commandline, but sent in as a list.
#' For example,
#' \Dontrun{
#' args <- check_args(arg_parser(c(
#'    "--config=gen_scaled_ar/rc_kappa.toml",
#'    "--country=gmb",
#'    "--years=2010:2011",
#'    "--overwrite"
#'    )))
#' args$country <- "uga"
#' funcmain(args)
#' }
#' @export
funcmain <- function(args) {
    # We will want to split this work different ways for development,
    # laptops, and the cluster, so we use an explicit domain decomposition.
    rampdata::initialize_workflow(args$config)
    options <- configr::read.config(args$config)[["options"]]

    # What can be done and which part we choose to do.
    available <- available_data(args$inversion, args$country, args$years)
    # How to split that up.
    plan <- plan_domain_decomposition(available, options)
    remaining_tiles <- remove_tiles_over_water(plan)

    # The part that this particular process will do.
    work <- task_work_from_tiles(remaining_tiles)
    flog.debug(paste("tile count", dim(work)[2], "size",
        paste(plan$tiles$blocksize, collapse = ",")))

    # By this point, someone told us what to do, so load all the data.
    process_extent <- raster_extent_from_work(work, plan$tiles$blocksize)
    flog.debug(paste(
        "extent", paste(process_extent, collapse = ","), "names",
        paste(names(process_extent), collapse = ","), "\n"))
    # The load_extent is what the process needs, relative to the domain coordinates.
    load_extent <- apply_relative_extent(plan$domain_extent, process_extent)
    data <- load_data(args$config, args$pr2ar, load_extent, args$years)
    # We permute data axes for more efficient running and
    # split the data into chunks for processing.
    timeseries_data <- prepare_timeseries(data, work, process_extent, plan$tiles$blocksize)
    data <- NULL  # Help the memory cleanup. This data is now in timeseries chunks.

    params <- timeseries_data$parameters
    # Uses the parallel random generation from the `parallel` package.
    set.seed(params$random_seed, "L'Ecuyer")
    draw_params <- draw_parameters(params, args$draws)
    chunks <- timeseries_data$chunks

    core_cnt <- parallel_core_cnt(args)
    flog.debug("%d chunks to %d cores", length(chunks), core_cnt)
    # Work in parallel on this machine.
    # Thread output usually goes to /dev/null, but we can redirect
    # it to a file. Put it with our data for now. Could be a mess
    # for thread contention, but they don't write a lot.
    thread_out_dir <- build_outvars_dir(args$outvars)
    thread_out <- rampdata::add_path(thread_out_dir, file = "threads.txt")
    thread_fn <- rampdata::as.path(thread_out)
    flog.debug(sprintf("thread file out to %s", thread_fn))
    if (core_cnt > 1) {
    cluster <- parallel::makeCluster(
    	    core_cnt,
	    type = "FORK",
	    outfile = thread_fn
	    )
    output <- parallel::clusterApplyLB(
        cluster,
        chunks,
        function(chunk) {
            flog.debug(".")
            over_block_draw(
                chunk,
                draw_params,
                timeseries_data$pr_to_ar_dt,
                params$confidence_percent
                )
        }
    )
    flog.debug("End of parallel work.")
    parallel::stopCluster(cluster)
    } else {
    output <- lapply(
        chunks,
        function(chunk) {
            flog.debug(".")
            over_block_draw(
                chunk,
                draw_params,
                timeseries_data$pr_to_ar_dt,
                params$confidence_percent
                )
        }
    )
    }

    # The output chunks need to be reassembled before writing.
    ready_to_write <- combine_output(output, plan$tiles$blocksize)
    write_output(ready_to_write, args$years, plan$domain_dimensions,
      plan$domain_extent, args, options)
}


#' If you want to run from the Bash command line, you can call this.
#' 
#' @export
main <- function() {
    flog.threshold(DEBUG)
    improved_errors()
    args <- check_args(arg_parser())
    funcmain(args)
}


#' Makes an initial plan.
#' @param args Command-line arguments. Doesn't care about the number of tasks,
#'     but otherwise might as well use the same command-line arguments here
#'     as elsewhere.
#' 
#' If you want to call this, see \code{\link{main}} for an example.
#' @export
construct_plan <- function(args) {
  # We will want to split this work different ways for development,
  # laptops, and the cluster, so we use an explicit domain decomposition.
  rampdata::initialize_workflow(args$config)
  options <- configr::read.config(args$config)[["options"]]

  # What can be done and which part we choose to do.
  available <- available_data(args$inversion, args$country, args$years)
  # How to split that up.
  plan <- plan_domain_decomposition(available, options)
  remaining_tiles <- remove_tiles_over_water(plan)

  # The part that this particular process will do.
  work <- task_work_from_tiles(remaining_tiles)
  flog.debug(paste("tile count", dim(work)[2], "size",
                   paste(plan$tiles$blocksize, collapse = ",")))
  out_dir <- build_outvars_dir(args$outvars)
  data.table::fwrite(t(work), file = rampdata::as.path(rampdata::add_path(out_dir, file = "work.csv")))
  save_plan(plan, rampdata::as.path(rampdata::add_path(out_dir, file = "plan.json")))
}


#' Reads the plan and produces this task's part of that plan.
#' @param args Command-line arguments. If you don't give it a task, it
#'     reads that task ID from the `SGE_TASK_ID` environment variable.
#' 
#' If you want to call this, see \code{\link{main}} for an example.
#' @export
worker <- function(args) {
  # We will want to split this work different ways for development,
  # laptops, and the cluster, so we use an explicit domain decomposition.
  rampdata::initialize_workflow(args$config)
  version_dir <- build_outvars_dir(args$outvars)
  plan_rp <- rampdata::add_path(version_dir, file = "plan.json")
  plan <- load_plan(rampdata::as.path(plan_rp))
  work_df <- data.table::fread(
    rampdata::as.path(rampdata::add_path(version_dir, file = "work.csv")))

  # The part that this particular process will do.
  my_work <- my_tasks(args$task, args$tasks, work_df)
  if (is.null(my_work)) {
    flog.info("Quitting because there is no work to do.")
    quit(save = "no", status = 0)
  }

  flog.debug(paste("tile count", dim(my_work)[2], "size",
                   paste(plan$tiles$blocksize, collapse = ",")))

  # By this point, someone told us what to do, so load all the data.
  process_extent <- raster_extent_from_work(my_work, plan$tiles$blocksize)
  flog.debug(paste(
    "extent", paste(process_extent, collapse = ","), "names",
    paste(names(process_extent), collapse = ","), "\n"))
  # The load_extent is what the process needs, relative to the domain coordinates.
  load_extent <- apply_relative_extent(plan$domain_extent, process_extent)
  data <- load_data(args$config, args$pr2ar, load_extent, args$years)
  # We permute data axes for more efficient running and
  # split the data into chunks for processing.
  timeseries_data <- prepare_timeseries(data, my_work, process_extent, plan$tiles$blocksize)
  data <- NULL  # Help the memory cleanup. This data is now in timeseries chunks.

  params <- timeseries_data$parameters
  # Uses the parallel random generation from the `parallel` package.
  set.seed(params$random_seed, "L'Ecuyer")
  draw_params <- draw_parameters(params, args$draws)
  chunks <- timeseries_data$chunks

  core_cnt <- parallel_core_cnt(args)
  flog.debug("%d chunks to %d cores", length(chunks), core_cnt)
  # Work in parallel on this machine.
  # Thread output usually goes to /dev/null, but we can redirect
  # it to a file. Put it with our data for now. Could be a mess
  # for thread contention, but they don't write a lot.
  thread_out_dir <- build_outvars_dir(args$outvars)
  thread_out <- rampdata::add_path(thread_out_dir, file = "threads.txt")
  thread_fn <- rampdata::as.path(thread_out)

  output <- lapply(
    chunks,
    function(chunk) {
      flog.debug(".")
      over_block_draw(
        chunk,
        draw_params,
        timeseries_data$pr_to_ar_dt,
        params$confidence_percent
      )
    }
  )
  chunks_rp <- rampdata::add_path(version_dir, file = sprintf("chunks%d.h5", args$task))
  chunks_fn <- rampdata::as.path(chunks_rp)
  save_outputs(chunks_fn, output)
}


#' Reads output from workers and produces GeoTIFFs and images.
#' @param args Command-line arguments. Ignores the task but uses the number of tasks.
#' 
#' If you want to call this, see \code{\link{main}} for an example.
#' @export
assemble <- function(args) {
  # We will want to split this work different ways for development,
  # laptops, and the cluster, so we use an explicit domain decomposition.
  rampdata::initialize_workflow(args$config)
  options <- configr::read.config(args$config)[["options"]]
  version_dir <- build_outvars_dir(args$outvars)
  plan_rp <- rampdata::add_path(version_dir, file = "plan.json")
  plan <- load_plan(rampdata::as.path(plan_rp))
  stopifnot(is_tile(plan$tiles$blocksize))

  task_name_fn <- function(task_idx) {
    chunks_rp <- rampdata::add_path(version_dir, file = sprintf("chunks%d.h5", task_idx))
    rampdata::as.path(chunks_rp)
  }

  missing <- list()
  for (miss_idx in 1:args$tasks) {
    find_fn <- task_name_fn(miss_idx)
    if (!file.exists(find_fn)) {
      missing[as.character(miss_idx)] <- find_fn
    }
  }
  if (length(missing) > 0) {
    flog.error(paste(
      sprintf("Cannot find the following %d tasks:", length(missing)),
      paste0(names(missing), collapse = ",")
      ))
    for (miss in missing) {
      flog.error(miss)
    }
    stopifnot(length(missing) == 0)
  }

  ds_names <- read_dataset_names(task_name_fn(1))
  # Subset to aeir and rc.
  subset <- NULL
  if (!is.null(subset) && length(subset) > 0) {
      keep_names <- vector(mode = "character", length = 0)
      for (var_keep in subset) {
          search_pattern <- sprintf("^%s", var_keep)
          keep_names <- c(keep_names, ds_names[grep(search_pattern, ds_names)])
      }
      best_names <- unique(keep_names)
  } else {
      best_names <- ds_names
  }

  core_cnt <- parallel_core_cnt(args)
  flog.info(sprintf("using %d cores", core_cnt))
  # To debug the parallel part, add outfile = "zrc.txt".
  cl <- parallel::makeCluster(core_cnt, outfile = "zrc.txt")
  doParallel::registerDoParallel(cl)
  # These may work better on Windows.
  # cl <- makeSOCKcluster(core_cnt)
  # registerDoSNOW(cl)

  foreach::foreach(
    loop_idx = seq_along(best_names),
    .packages = c("globalrc")
    ) %dopar% {
      ds_name <- best_names[loop_idx]
      output <- list()
      # Rampdata keeps the config in package namespace, so that has to be
      # initialized in the workers.
      rampdata::initialize_workflow(args$config)
      for (task_idx in 1:args$tasks) {
        more_output <- read_outputs(task_name_fn(task_idx), ds_name)
        output <- c(output, more_output)
      }
      # The output chunks need to be reassembled before writing.
      ready_to_write <- combine_output(output, plan$tiles$blocksize, ds_name)
      write_output(
        ready_to_write, args$years, plan$domain_dimensions, plan$domain_extent, args, options
        )
  }
  
  parallel::stopCluster(cl)
}
