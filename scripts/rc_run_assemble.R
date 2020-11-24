library(globalrc)

flog.threshold(DEBUG)
globalrc::improved_errors()
args <- globalrc::check_args(globalrc::arg_parser())
globalrc::assemble(args)
