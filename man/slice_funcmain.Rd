% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slice_main.R
\name{slice_funcmain}
\alias{slice_funcmain}
\title{If you want to run from the R command line, you can call this.}
\usage{
slice_funcmain(args)
}
\arguments{
\item{args}{Command-line arguments, after they've been checked.

Use the same arguments as the commandline, but sent in as a list.
For example,\preformatted{args <- check_args(arg_parser(c(
   "--config=gen_scaled_ar/rc_kappa.toml",
   "--country=gmb",
   "--years=2010:2011",
   "--overwrite"
   )))
args$country <- "uga"
funcmain(args)
}}
}
\description{
If you want to run from the R command line, you can call this.
}
