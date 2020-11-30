
#' Check that this is a tile.
#' @param p A vector that should have names `row` and `col`
#' @return A bool for yes or no.
#' These are the two main data structures. Tiles and blocksizes are measured
#' in rows and columns, and they are stored as a numeric vector with names
#' "row" and "col". If you forget those, then you get NA. The same goes for
#' extents, which all have "rmin", "rmax", "cmin", "cmax" for the mininum and
#' maximum of rows and columns. We use these names to make it easier to
#' visually ascertain that rows are matched with rows and columns with columns.
is_tile <- function(p) all(names(p) == c("row", "col"))
is_extent <- function(e) all(names(e) == c("rmin", "rmax", "cmin", "cmax"))


#' Compute a domain extent from a process extent within a domain extent.
#' @param domain_extent The larger bounding box.
#' @param process_extent A bounding box within the domain extent.
#' @return A bounding box measured from the domain extent's coordinate system.
apply_relative_extent <- function(domain_extent, process_extent) {
    stopifnot(is_extent(domain_extent))
    stopifnot(is_extent(process_extent))
    offset <- domain_extent[c("rmin", "cmin")] - 1
    process_extent + unname(offset[c("rmin", "rmin", "cmin", "cmin")])
}


#' Find where the tile is within the process extent.
#' @param process_extent An extent vector representing what tiles this process
#'     neeeds to load.
#' @param tile_extent An extent vector which is in numbers of tiles.
#' @return Where the tiles are within the process.
#' Both are relative to the domain extent.
find_relative_extent <- function(process_extent, tile_extent) {
    stopifnot(is_extent(process_extent))
    stopifnot(is_extent(tile_extent))
    offset <- process_extent[c("rmin", "cmin")] - 1
    tile_extent - unname(offset[c("rmin", "rmin", "cmin", "cmin")])
}


#' Removes tiles that have all NA for this example slice.
#' @param available This is a description of where all tiles are.
#' @return A list of remaining tiles, by row and column, so a 2D array.
#' Tile is within the domain. The sample data is already cut
#' to that domain.
remove_tiles_over_water <- function(available) {
    stopifnot("sample_slice" %in% names(available))
    # This slice is cut from the domain.
    pfpr <- available$sample_slice
    tile_cnt <- available$tiles$tile_cnt
    blocksize <- available$tiles$blocksize

    ridx <- rep(1:tile_cnt["row"], tile_cnt["col"])
    cidx <- rep(1:tile_cnt["col"], each = tile_cnt["row"])
    flog.debug("examining %d tiles for NA", length(ridx))
    tiles <- rbind(ridx, cidx)
    # If there is only one tile, we need to insist this is 2D array.
    dim(tiles) <- c(2, length(ridx))
    dimnames(tiles) <- list(c("row", "col"), NULL)
    keep_tile <- vapply(1:(dim(tiles)[2]), function(tile_idx) {
        extent <- tile_extent(tiles[, tile_idx], blocksize)
        rmax <- min(extent["rmax"], dim(pfpr)[1])
        cmax <- min(extent["cmax"], dim(pfpr)[2])
        chunk <- pfpr[
            extent["rmin"]:rmax,
            extent["cmin"]:cmax
        ]
        sum(!is.na(chunk)) > 0
    }, FUN.VALUE = logical(1))
    removed_cnt <- sum(!keep_tile)
    kept_cnt <- sum(keep_tile)
    flog.info(sprintf("Removing %d tiles with no data keep %d",
        removed_cnt, kept_cnt))
    remaining <- tiles[, keep_tile]
    # Again, if there is only one to keep, insist it be 2D.
    dim(remaining) <- c(2, kept_cnt)
    dimnames(remaining) <- list(c("row", "col"), NULL)
    remaining
}


tile_extent <- function(tile, blocksize) {
    stopifnot(is_tile(tile))
    stopifnot(is_tile(blocksize))
    c(
        rmin = unname((tile["row"] - 1L) * blocksize["row"] + 1L),
        rmax = unname(tile["row"] * blocksize["row"]),
        cmin = unname((tile["col"] - 1L) * blocksize["col"] + 1L),
        cmax = unname(tile["col"] * blocksize["col"])
    )
}


#' Given a list of tiles, figures out the bounding box of pixel values.
#' @param work A 2d matrix with a dimname on the first dimension of "row", "col".
#' @param blocksize A row and col size for blocks for the work.
raster_extent_from_work <- function(work, blocksize) {
    stopifnot(all(dimnames(work)[[1]] == c("row", "col")))
    stopifnot(is_tile(blocksize))
    small <- tile_extent(c(row = min(work["row",]), col = min(work["col", ])), blocksize)
    large <- tile_extent(c(row = max(work["row",]), col = max(work["col", ])), blocksize)
    c(small["rmin"], large["rmax"], small["cmin"], large["cmax"])
}
