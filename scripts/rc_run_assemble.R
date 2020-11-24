#pkgload::load_all("..")
# From r-lib/testthat in the test-files.R file.
# pkgload::load_all(test_dir, helpers = FALSE, quiet = TRUE)
library(futile.logger)
library(globalrc)

flog.threshold(DEBUG)
globalrc::improved_errors()
args <- globalrc::check_args(globalrc::arg_parser())
globalrc::assemble(args)
