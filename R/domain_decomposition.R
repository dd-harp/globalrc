
#' Given command-line arguments, determine how much data to load.
#' It's deciding splits without consulting arguments for now.
#' This decomposes the domain into tiles because, when we parallelize
#' with the parallel package, each task needs at least a few seconds
#' of work to do in order to make data transfer worthwhile.
#' @param available A description of what data could be processed.
#' @param options These parameters define `blocksize` and `single_tile_max`.
#' @return The same list with an added `tiles` member that describes
#'     the dimensions of computational tiles relative to the origin
#'     of the `domain_extent`, the portion we will compute.
plan_domain_decomposition <- function(available, options) {
    extent <- available$domain_extent
    stopifnot(is_extent(extent))
    rows <- (extent["rmax"] - extent["rmin"] + 1L)
    cols <- (extent["cmax"] - extent["cmin"] + 1L)
    if (rows * cols < options$single_tile_max) {
        tile_cnt <- c(row = 1L, col = 1L)
        blocksize <- c(row = unname(rows), col = unname(cols))
    } else {
        bs <- as.integer(options$blocksize)
        blocksize <- c(row = bs, col = bs)
        tile_cnt <- c(
            row = as.integer(ceiling(rows / blocksize["row"])),
            col = as.integer(ceiling(cols / blocksize["col"]))
            )
    }
    available[["tiles"]] <- list(
        tile_cnt = tile_cnt,
        blocksize = blocksize
    )
    available
}


#' Takes the extent of the tiles and returns a list of which tiles to do.
#' @param tiles A list of tiles as a 2xN matrix.
#' @return An array that's 2 x the number of tiles. It's row, col.
#' Assumes that all work is done by the same process.
#' This orders the work so that it uses cache better.
task_work_from_tiles <- function(tiles) {
    # The Hilbert curve is a space-filling curve, so it orders the 2D tiles.
    tidx <- vapply(1:dim(tiles)[2], FUN = function(tile_idx) {
            t <- tiles[, tile_idx]
            encode_hilbert(t["row"], t["col"])
        }, FUN.VALUE = integer(1)
        )
    by_curve <- order(tidx)
    tiles[, by_curve]
    dim(tiles) <- c(2, length(by_curve))
    dimnames(tiles) <- list(c("row", "col"), NULL)
    tiles
}


#' Takes the extent of the tiles and returns a list of which tiles to do.
#' @param tile_cnt The tile_cnt from the plan from `plan_domain_decomposition`.
#' @return An array that's 2 x the number of tiles. It's row, col.
#' Assumes that all work is done by the same process.
#' This orders the work so that it uses cache better.
task_work <- function(tile_cnt) {
    stopifnot(is_tile(tile_cnt))
    ridx <- rep(1:tile_cnt["row"], tile_cnt["col"])
    cidx <- rep(1:tile_cnt["col"], each = tile_cnt["row"])
    # The Hilbert curve is a space-filling curve, so it orders the 2D tiles.
    tidx <- vapply(1:length(cidx), FUN = function(tile_idx) {
            encode_hilbert(ridx[tile_idx], cidx[tile_idx])
        }, FUN.VALUE = integer(1)
        )
    by_curve <- order(tidx)
    work <- rbind(ridx[by_curve], cidx[by_curve])
    dimnames(work) <- list(c("row", "col"), NULL)
    work
}


#' Pull my work, according to my task index.
#' @param task_id integer index of task within all tasks. One-based.
#' @param task_cnt integer total number of tasks.
#' @param work Data frame with row and column of tile.
#' @return 2D array that is (row, col) by work items.
#'
#' The math here is so that, if the tasks don't split evenly, there
#' isn't a single task with only one thing to do. The first few
#' tasks get one more, and the rest get one less.
my_tasks <- function(task_id, task_cnt, work) {
  n <- nrow(work)
  m <- task_cnt
  w <- 1 + (n - 1) %/% m
  x <- n - m * (w - 1)
  if (task_id <= x) {
    rows <- (1 + (task_id - 1) * w):(task_id * w)
  } else {
    if (w > 1) {
      r <- task_id - x
      b <- x * w
      rows <- (b + 1 + (r - 1) * (w - 1)):(b + r * (w - 1))
    } else {
      return(NULL)
    }
  }
  as.matrix(t(work[rows, ]))
}
