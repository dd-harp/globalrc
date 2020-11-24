
#' Parse command-line arguments.
#' @param args Optionally pass the results of `commandArgs(TRUE)`.
#'     We have this parameter so that we can test without side effects.
#' @return A list where args that aren't on the command line are NULL.
#' @export
arg_parser <- function(args = NULL) {
    doc <- "pr to Rc

Usage:
  rc_kappa.R [options]
  rc_kappa.R (-h | --help)

Options:
  -h --help              Show help.
  --config=<config>      A configuration file.
  --country=<alpha3>     The three-letter country code for a country's outline.
  --outvars=<outversion> Version of output to write.
  --overwrite            Whether to overwrite outputs.
  --years=<year_range>   A range of years to do, as an R range, 2000:2010.
  --cores=<core_cnt>     Tell it how many cores to use in parallel.
  --draws=<draw_cnt>     How many draws to use.
  --task=<task_id>       If this is a task, which task.
  --tasks=<task_cnt>     Total number of tasks. You should set this for workers.
"
    if (is.null(args)) {
        args <- commandArgs(TRUE)
    }
    parsed_args <- docopt::docopt(doc, version = "rc_kappa 1.0", args = args)
    if (is.null(parsed_args$config)) {
        parsed_args$config <- "rc_kappa.toml"
    }
    parsed_args
}


#' Check whether the input arguments make sense.
#' @param args A list of command-line arguments resulting from the arg parser.
#' @return The same list of command-line arguments with some new types.
#' @export
check_args <- function(args) {
    OK <- TRUE
    if (!file.exists(args$config)) {
        flog.error(paste("Cannot find file", args$config, "from the config argument.",
            "from the directory", getwd()))
        OK <- FALSE
    } else {
        if (!is.list(configr::read.config(args$config))) {
            flog.error(paste("configr says the config file has a syntax error", args$config))
            OK <- FALSE
        }
    }
    if (!is.null(args$country) && nchar(args$country) != 3) {
        flog.error(paste("The country should be the three-letter Alpha-3 code:", args$country))
        OK <- FALSE
    }
    if (!is.null(args$years)) {
        start_year <- suppressWarnings(as.integer(substr(args$years, 1, 4)))
        end_year <- suppressWarnings(as.integer(substr(args$years, 6, 9)))
        if (is.na(start_year) | is.na(end_year)) {
            flog.error(paste("Expecting the year range to be 2000:2010. Found", args$years))
            OK <- FALSE
        } else {
            args$years <- start_year:end_year
        }
    }
    # The task ID is usually not set. If it is, then we take it.
    # If it isn't, we look at the SGE_TASK_ID environment variable.
    if (!is.null(args$task)) {
      args$task <- tryCatch(as.integer(args$task), warning = function(x) {
        stop(paste("Could not convert task id to integer", args$task))
      })
    } else {
      task_id <- Sys.getenv("SGE_TASK_ID", unset = NA)
      if (!is.na(task_id)) {
        args$task <- as.integer(task_id)
      }
    }
    if (!is.null(args$tasks)) {
      args$tasks <- tryCatch(as.integer(args$tasks), warning = function(x) {
        stop(paste("Could not convert tasks cnt to integer", args$tasks))
      })
    }
    args$cores <- parallel_core_cnt(args)
    stopifnot(OK)
    args
}
