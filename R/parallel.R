parallel_core_cnt <- function(args = NULL) {
    # So the user can specify in command-line arguments.
    if (is.list(args) && "cores" %in% names(args) && !is.null(args$cores)) {
        cores <- suppressWarnings(as.integer(args$cores))
        if (!is.na(cores) && length(cores) == 1L && cores > 0L) {
            return(cores)
        } else {
            flog.warn(paste("Could not use args$cores to set core count:", args$cores))
        }
    }
    # Because we run on the cluster, and the cluster sets environment variables.
    for (env_var in c("NCPUS", "SGE_HGR_fthread")) {
        env_str <- Sys.getenv(env_var, unset = NA)
        if (!is.na(env_str)) {
            env_cores <- suppressWarnings(as.integer(env_str))
            if (length(env_cores) == 1L && env_cores > 0L) {
                return(env_cores)
            }
        }
    }
    # Because we want all the cores otherwise.
    core_cnt <- parallel::detectCores()
    if (is.integer(core_cnt)) {
        return(core_cnt)
    } else {
        return(2L)
    }
}
