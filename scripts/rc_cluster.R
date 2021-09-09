# This script builds shell scripts to submit to the cluster
# in order to do a run of globalrc.

# This libpath is because I'm using renv during development, so globalrc
# is in this little local library, and it isn't found when I run R
# from the scripts directory.
has_globalrc <- requireNamespace("globalrc", quietly=TRUE)
if (!has_globalrc) {
  .libPaths(c(.libPaths(), "../renv/library/R-4.0/x86_64-pc-linux-gnu"))
}

library(futile.logger)
library(docopt)
library(globalrc)

cluster_arg_parser <- function(args = NULL) {
  doc <- "pr to Rc

Usage:
  rc_cluster.R [options]
  rc_cluster.R (-h | --help)

Options:
  -h --help              Show help.
  --config=<config>      A configuration file to use as input.
  --outvars=<outversion> Version of output to write.
  --years=<year_range>   A range of years to do, as an R range, 2000:2010.
  --draws=<draw_cnt>     How many draws to use.
  --country=<alpha3>     The three-letter country code for a country's outline.
  --tasks=<task_cnt>     Total number of tasks. You should set this for workers.

You need to give it the config, outvars, years, and draws.
The country and tasks are optional. If you don't specify tasks, it will
request a number of tasks to put in the worker.
"
  if (is.null(args)) {
    args <- commandArgs(TRUE)
  }
  parsed_args <- docopt::docopt(doc, version = "rc_kappa 1.0", args = args)
  parsed_args$config <- normalizePath(parsed_args$config, mustWork = FALSE)
  parsed_args
}


# Examine the data to find size of tiles.
make_plan <- function() {
  flog.info("make_plan")
  plan_args <- globalrc::check_args(globalrc::arg_parser())
  tile_cnt <- construct_plan(plan_args)
  tile_cnt
}


# These are values to replace in template scripts.
build_replaces <- function(args, task_cnt) {
  flog.info("build_replaces")
  if (!is.null(args$country)) {
    country_str <- sprintf("--country=%s", args$country)
  } else {
    country_str <- ""
  }

  list(
    tasks = task_cnt,
    config = args$config,
    outvars = args$outvars,
    years = args$years,
    draws = args$draws,
    country = country_str,
    MISSING = "${MISSING}",  # These replace with the same values.
    SGE_TASK_ID = "${SGE_TASK_ID}"
  )
}


# If the person specified it, then don't go get one.
get_task_cnt <- function(tile_cnt, task_cnt) {
  flog.info("get_task_cnt")
  if (is.null(task_cnt)) {
    msg <- sprintf("How many tasks for %d tiles? ", tile_cnt)
    #task_cnt_str <- readline(prompt = msg)
    cat(msg)
    task_cnt_str <- readLines("stdin", n = 1)
    task_cnt <- as.integer(task_cnt_str)
  }  # else we were given the task count
  task_cnt
}


# Reads a template and writes something with the version name.
make_script <- function(replaces, role) {
  flog.info("make_script")
  base <- sprintf("rc_%s_base.sh", role)
  base_worker <- readLines(base)
  write_it <- with(replaces, {
    outlines <- vector(mode = "character", length = length(base_worker))
    for (lidx in seq_along(base_worker)) {
      outlines[lidx] <- stringr::str_interp(base_worker[lidx])
    }
    outlines
  })
  out_fn <- sprintf("%s_%s.sh", replaces$outvars, role)
  flog.info(sprintf("Writing %s to %s", role, out_fn))
  writeLines(write_it, out_fn)
  out_fn
}


main <- function() {
  flog.threshold(INFO)
  args <- cluster_arg_parser()
  tile_cnt <- make_plan()
  task_cnt <- get_task_cnt(tile_cnt, args$tasks)
  replaces <- build_replaces(args, task_cnt)
  worker <- make_script(replaces, "worker")
  assemble <- make_script(replaces, "assemble")
  cat("Run with:\n")
  cat(sprintf("qsub %s\n", worker))
  cat(sprintf("qsub %s -hold_jid <workerid>\n", assemble))
}

main()
