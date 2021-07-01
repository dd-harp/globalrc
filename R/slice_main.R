# This is a main script that processes maps one time-slice at a time.
# This should be simpler to modify than other parallelizing strategies.

library(futile.logger)


#' If you want to run from the R command line, you can call this.
#' @param args Command-line arguments, after they've been checked.
#'
#' Use the same arguments as the commandline, but sent in as a list.
#' For example,
#'
#' ```
#' args <- check_args(arg_parser(c(
#'    "--config=gen_scaled_ar/rc_kappa.toml",
#'    "--country=gmb",
#'    "--years=2010:2011",
#'    "--overwrite"
#'    )))
#' args$country <- "uga"
#' funcmain(args)
#' ```
#'
#' @export
slice_funcmain <- function(args) {
  rampdata::initialize_workflow(args$config)
  options <- configr::read.config(args$config)[["options"]]

  # What can be done and which part we choose to do.
  available <- available_data(args$inversion, args$country, args$years)
  load_extent <- available$domain_extent
  years <- available$year[1]
  data <- load_data(args$config, args$pr2ar, load_extent, years)

  params <- data$parameters
  # Uses the parallel random generation from the `parallel` package.
  set.seed(params$random_seed, "L'Ecuyer")
  draws <- ifelse(is.null(args$draws), 100, args$draws)
  draw_params <- draw_parameters(params, draws)
  thread_out_dir <- globalrc:::build_outvars_dir(args$outvars)
  chunk <- list(
    am = data$am$`2000`,
    pfpr = data$pfpr$`2000`
  )
  list_of_summaries <- over_block_draw(
    chunk,
    draw_params,
    data$pr_to_ar_dt,
    params$confidence_percent
  )
  multiyear <- lapply(list_of_summaries, function(s) {
    dim(s) <- c(1, dim(s))
    s
  })
  write_output(
    multiyear, years, available$domain_dimensions, load_extent,
    args, options
    )
}


#' If you want to run from the Bash command line, you can call this.
#'
#' @export
slice_main <- function() {
  flog.threshold(DEBUG)
  improved_errors()
  args <- check_args(arg_parser())
  slice_funcmain(args)
}
