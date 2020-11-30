#' These functions are used to arrange a walk across a 2D matrix to 
#' use cache more efficiently.
#' They implement a space-filling curve in 2D.
#' A space-filling curve gives you a single integer to index a two-dimensional
#' array. If you have 2D tiles, then you would find the Hilbert
#' index of each tile and sort by that index in order to improve cache use.


#' Given an x and y value, translate it into a Hilbert z-value.
#' @param x An integer for one axis.
#' @param y An integer for the other axis.
#' @return A z-value to index the two of them.
#' @export
encode_hilbert <- function(x, y) {
    encode_hilbert_zero(x - 1L, y - 1L) + 1L
}


#' Given a Hilbert z-value, translate it into an x and y value.
#' @param z An integer to index a 2D matrix element.
#' @return An array `c(x, y)` of array coordinates.
#' @export
decode_hilbert <- function(z) {
    hilbert <- decode_hilbert_zero(z - 1L)
    c(hilbert[1] + 1L, hilbert[2] + 1L)
}


largest_hilbert_z <- function(x, y) {
    zmax <- 0
    j <- x
    for (i in 1:y) {
        zmax <- max(zmax, encode_hilbert(i, j))
    }
    i <- y
    for (j in 1:x) {
        zmax <- max(zmax, encode_hilbert(i, j))
    }
    zmax
}


.hilbert.variant.encode.table = matrix(
    c(
    function(x, y, w) { c(x, y) },           function(x, y, w) { c(x, y) },
    function(x, y, w) { c(y-w, x) },         function(x, y, w) { c(y, x-w) },
    function(x, y, w) { c(y-w, x-w) },       function(x, y, w) { c(y-w, x-w)},
    function(x, y, w) { c(2L*w-x-1L, w-y-1L)}, function(x, y, w) {c(w-x-1L, 2L*w-y-1L)}),
    nrow = 4,
    ncol = 2,
    byrow = TRUE
)

#' Computes an integer Hilbert index for x and y using a variant algorithm.
#' 
#' @param x An integer x-axis value, 0-based.
#' @param y An integer y-axis value, 0-based.
#' @return A hilbert z-value, an integer.
#' 
#' Given two integer indices for a 2-dimensional plane, return a single index.
#' This index is designed to increase locality for 1-dimensional access.
#' It does this by keeping nearby points in 2 dimensions also nearby in
#' 1 dimension.
#' The variant algorithm used differs from the usual Hilbert code because it
#' doesn't need to know the size of the whole grid before computing the code.
#' It looks like a slightly-rotated version of the Hilbert curve, but it
#' has the benefit that it is 1-1 between `(x, y)` and `z`, so you can translate
#' back and forth.
#' This function is zero-based. `0 <= x < 2^n`, `0 <= y < 2^n`, and the result
#' is `0 <= z < 4^n`.
#'
#' N. Chen, N. Wang, B. Shi, A new algorithm for encoding and decoding the
#' Hilbert order. Software—Practice and Experience 2007; 37(8): 897–908.
#'
#' See also: [`decode_hilbert_zero`](@ref), [`encode_hilbert`](@ref).
#' @export
encode_hilbert_zero <- function(x, y) {
    z <- 0L
    if (x == 0L & y == 0L) {
        return(z)
    }
    rmin <- as.integer(floor(log2(max(x, y))) + 1L)
    w <- bitwShiftL(1, rmin - 1L)
    while (rmin > 0L) {
        if (bitwAnd(rmin, 1L) == 1L) {
            quadrant <- ifelse(x < w, ifelse(y < w, 0L, 1L), ifelse(y < w, 3L, 2L))
            parity <- 1L
        } else {
            quadrant <- ifelse(x < w, ifelse(y < w, 0L, 3L), ifelse(y < w, 1L, 2L))
            parity <- 2L
        }
        z <- bitwShiftL(z, 2L) + quadrant
        xy <- .hilbert.variant.encode.table[quadrant + 1, parity][[1]](x, y, w)
        x <- xy[1]
        y <- xy[2]
        rmin <- rmin - 1L
        w <- bitwShiftR(w, 1L)
    }
    z
}


.hilbert.variant.decode.table = matrix(c(
    function(x, y, w) {c(x, y)},    function(x, y, w) {c(x, y)},
    function(x, y, w) {c(y, x+w)},  function(x, y, w) {c(y+w, x)},
    function(x, y, w) {c(y+w, x+w)}, function(x, y, w) {c(y+w, x+w)},
    function(x, y, w) {c(2*w-x-1, w-y-1)}, function(x, y, w) {c(w-x-1, 2*w-y-1)}),
    nrow = 4,
    ncol = 2,
    byrow = TRUE
)


decode_hilbert_zero <- function(z) {
    r <- bitwAnd(z, 3L)
    xy_init <- list(c(0L, 0L), c(0L, 1L), c(1L, 1L), c(1L, 0L))[[r + 1L]]
    x <- xy_init[1]
    y <- xy_init[2]
    z <- bitwShiftR(z, 2L)
    rmin <- 2L
    w <- 2L
    while (z > 0L) {
        r <- bitwAnd(z, 3L)
        parity <- 2L - bitwAnd(rmin, 1L)
        xy <- .hilbert.variant.decode.table[r + 1L, parity][[1]](x, y, w)
        x <- xy[1]
        y <- xy[2]
        z <- bitwShiftR(z, 2L)
        rmin <- rmin + 1L
        w <- bitwShiftL(w, 1L)
    }
    c(x, y)
}


test_hilbert <- function() {
    for (i in 1:400) {
        z <- as.integer(i - 1L)
        xy <- decode_hilbert_zero(z)
        zprime <- encode_hilbert_zero(xy[1], xy[2])
        stopifnot(z == zprime)
    }
}
