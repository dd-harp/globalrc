#' Find nonzeros in each chunk of tile data to see they loaded.
#' @param chunks A list of tiles with data.
chunk_nonzero <- function(chunks) {
    vapply(
        chunks,
        function(chunk) {
            c(sum(!is.na(chunk$pfpr)), sum(!is.na(chunk$am)))
        },
        integer(2)
    )
}


#' Find the sum of each chunk of tile data to see they differ.
#' @param chunks A list of tiles with data.
chunk_sums <- function(chunks) {
    vapply(
        chunks,
        function(chunk) {
            c(sum(chunk$pfpr, na.rm = TRUE), sum(chunk$am, na.rm = TRUE))
        },
        numeric(2)
    )
}



#' Check that the results are all finite.
#' @param results A list of arrays.
#' @param inputs The input list of arrays, used for debugging failure.
#' @return The input results. This will stop on failure.
#' We shouldn't have any NANs here. Good data in, good data out.
#' This is what let's us ignore NA later. This function is expensive
#' because it scans all points, but it's important.
check_linearized_outputs <- function(results, inputs) {
    for (check_idx in seq_along(results)) {
        check_arr <- results[[check_idx]]
        if (!all(is.finite(check_arr))) {
            bad_ones <- !is.finite(check_arr)
            res_bad <- data.frame(c(inputs, results))[bad_ones,]
            debug_file <- "debug.csv"
            data.table::fwrite(res_bad, file = debug_file)
            inf_cnt <- sum(bad_ones)
            nan_cnt <- sum(is.na(check_arr))
            var_name <- names(results)[check_idx]
            flog.warn(sprintf(
                "Found %d nan %d inf in %s. Writing to %s\n",
                nan_cnt, inf_cnt, var_name,
                normalizePath(debug_file)
                ))
            stopifnot(nan_cnt == 0 && inf_cnt == 0)
        }
    }
    results
}


#' Takes an input array in multiple dimensions and linearizes and removes NA.
#' @param input_list A list of arrays, one for each variable.
#' @param run_func The function to run over this list.
#' @return A list of arrays, the output of the `run_func`.
#' This replaces the outputs back into multiple dimensions after making the call.
linearized_work <- function(input_list, run_func) {
    # There are NA values in the arrays in the input_list.
    # The arrays are probably three-dimensional.
    not_available <- lapply(input_list, function(check) is.na(check))
    input_dims <- dim(input_list[[1]])
    stopifnot(length(input_dims) == 3)
    dim_keep <- length(input_dims)
    # If any array has an NA in a position in the array, then that position is
    # thrown out for all input arrays.
    keep <- !rowSums(array(do.call(c, not_available), dim = c(input_dims, 2)), dims = dim_keep)
    linear <- lapply(input_list, function(data) {
        as.numeric(data[keep])
    })

    # How could we have a tile with no data, when they were excluded?
    # They were excluded based on PR, not AM.
    if (sum(keep) > 0) {
        results <- check_linearized_outputs(run_func(linear), linear)
    } else {
        # Don't call the function with no data, but how do we know
        # variable names? Call the function with small, fake data.
        just_names <- run_func(list(pfpr=c(0.1), am = c(0.02)))
    results <- lapply(just_names, function(x) numeric(0))
    }

    # The returned values will have NA for the same pattern in all outputs.
    lapply(results, function(result) {
        complete <- array(dim = input_dims)
        complete[keep] <- result
        array(complete, dim = input_dims)
    })
}


#' The draws are hierarchical by variable and then quantile. This flattens them.
#' @param draw_list The hierarchical list with names.
#' @return A single list with the two names concatenated.
flatten_draws <- function(draw_list) {
    ll <- draw_list
    sl_cnt <- length(ll[[1]])
    draw_names <- lapply(1:(sl_cnt * length(ll)), function(idx) {
        list_idx <- 1 + (idx - 1) %/% sl_cnt
        sub_list <- ll[[list_idx]]
        sublist_idx <- 1 + (idx - 1) %% sl_cnt
        name <- names(ll)[list_idx]
        sub_name <- names(sub_list)[sublist_idx]
        paste(name, sub_name, sep = "_")
    })
    draw_vals <- lapply(1:(sl_cnt * length(ll)), function(idx) {
        list_idx <- 1 + (idx - 1) %/% sl_cnt
        sub_list <- ll[[list_idx]]
        sublist_idx <- 1 + (idx - 1) %% sl_cnt
        sub_list[[sublist_idx]]
    })
    names(draw_vals) <- draw_names
    draw_vals
}


#' Given a list of draws, summarize mean and confidence intervals.
#' @param draws A list of draws. Each draw is a list of variables,
#'     where each variable is a three-dimensional array.
#' @param confidence_percent As in, 95.
#' @return a list with one name for each variable. Each list
#'     item has three parts, lower, median, and upper, for that variable.
#'     Each of lower, median, and upper has dimensions year, row, col.
summarize_draws <- function(draws, confidence_percent) {
    ci <- 0.5 * (1 - 0.01 * confidence_percent)
    quantiles <- c(ci, 0.5, 1 - ci)
    draw_cnt <- length(draws)
    array_names <- names(draws[[1]])
    array_dim <- dim(draws[[1]][[1]])
    stopifnot(length(array_dim) == 3)
    summarized <- lapply(array_names, function(name) {
        var_data <- array(
            do.call(
                c,
                lapply(1:draw_cnt, function(combine_idx) {
                    draws[[combine_idx]][[name]]
                })
                ),
            dim = c(array_dim, draw_cnt)
        )
        draws_first <- aperm(var_data, c(4, 1, 2, 3))
        quantile_first <-apply(
            draws_first,
            MARGIN = c(2, 3, 4),
            # We can ignore NA here b/c we checked earlier, and some pixels
            # _should be all NA_ in water.
            FUN = function(x) stats::quantile(x, quantiles, na.rm = TRUE)
            )
        quantile_last <- aperm(quantile_first, c(2, 3, 4, 1))
        sumd <- list(
            lower = quantile_last[, , , 1, drop= FALSE],
            median = quantile_last[, , , 2, drop= FALSE],
            upper = quantile_last[, , , 3, drop= FALSE]
        )
        sumd <- lapply(sumd, function(x) {
            pre_dim <- dim(x)
            dim(x) <- pre_dim[1:3]
            x
        })
        stopifnot(length(dim(sumd[["median"]])) == 3)
        sumd
    })
    names(summarized) <- array_names
    flattened <- flatten_draws(summarized)
    lapply(flattened, function(x) {
      dimnames(x) <- list(year=NULL, row=NULL, col=NULL)
      x
    })
}


#' Run the pixel function on a chunk of data.
#' @param chunk a list with the pfpr and am datasets.
#' @param parameters a data.frame, at this point, with a row per draw.
#' @param pr_to_ar_dt The pr2ar mesh data.
#' @param confidence_percent Confidence interval, out of 100.
over_block_draw <- function(chunk, parameters, pr_to_ar_dt, confidence_percent) {
    stopifnot(is.numeric(confidence_percent))
    strategies <- list(
        pr_to_ar = ar_of_pr_rho(pr_to_ar_dt),
        ar2pr = build_ar2pr(pr_to_ar_dt)
    )
    # The function we call assumes that every entry in the chunk
    # is a dataset, so remove the tile.
    only_data <- chunk[names(chunk)[!names(chunk) %in% "tile"]]
    # The first dimension is the timeseries, so row=2, col=3.
    draw_cnt <- nrow(parameters)
    draws <- lapply(1:draw_cnt, function(draw_idx) {
        run_func <- function(plaquette) {
            pixel_four(plaquette$pfpr, plaquette$am, parameters[draw_idx, ], strategies)
        }
        linearized_work(only_data, run_func)
    })
    summarized <- summarize_draws(draws, confidence_percent)
    flog.debug(paste("summarized names", paste0(names(summarized), collapse = ",")))
    summarized[["block"]] <- chunk$tile
    summarized
}
