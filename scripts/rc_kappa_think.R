
pfpr_dir <- rampdata::ramp_path("/globalrc/inputs/PfPR_medians_2000_2019/201029")
am_dir <- rampdata::ramp_path("/globalrc/inputs/AM_1991_2002/201029")

pfpr_file <- rampdata::add_path(pfpr_dir, file = "PfPR_median_Africa_admin0_2000.tif")
file.exists(rampdata::as.path(pfpr_file))
pfpr_years <- years_in_filenames(list.files(rampdata::as.path(pfpr_dir), "*.tif"))

am_file <- rampdata::add_path(am_dir, file = "2000.effective.treatment.tif")
file.exists(rampdata::as.path(am_file))
am_years <- years_in_filenames(list.files(rampdata::as.path(am_dir), "*.tif"))

both_years <- sort(intersect(am_years, pfpr_years))

pfpr <- raster::raster(rampdata::as.path(pfpr_file))
am <- raster::raster(rampdata::as.path(am_file))
#plot(am)
am_mat <- raster::as.matrix(am)
pr_mat <- raster::as.matrix(pfpr)
stopifnot(all(dim(pr_mat) == dim(am_mat)))

library(raster)
live <- 1100000
prval <- getValues(pfpr)
stopifnot(prval[live] > 0)  # 0.7186912
row <- rowFromCell(pfpr, live)  # 655
col <- colFromCell(pfpr, live)  # 626
blocked <- raster::getValuesBlock(pfpr, row = 655, nrows = 1, col = 626, ncols = 1)
stopifnot(blocked == prval[live])
extracted <- raster::extract(pfpr, matrix(c(626, 655), ncol = 2))
stopifnot(extracted == prval[live])
stopifnot(pr_mat[655, 626] == prval[live])
# Now it gets surprising. The linear getValues isn't arranged the
# way you arrange R matrices. It's C order, not Fortran order.
stopifnot(prval[626 + (655 - 1) * raster::ncol(pfpr)] == prval[live])


r1000 <- matrix(
    raster::getValues(pfpr, row = 655, nrows = 2),
    nrow = 2,
    byrow = TRUE
)
stopifnot(dim(pfpr)[2] * 2 == length(r1000))
stopifnot(r1000[1, 626] > 0)
stopifnot(r1000[1, 626] == prval[live])

rblock1 <- matrix(
    raster::getValuesBlock(pfpr, row = 655, nrows = 2, col = 1, ncols = dim(pfpr)[2]),
    nrow = 2,
    byrow = TRUE
)
rblock <- matrix(
    raster::getValuesBlock(pfpr, row = 655, nrows = 2, col = 626, ncols = 3),
    nrow = 2,
    byrow = TRUE
)
stopifnot(all(rblock == rblock1[, 626:628]))


bounds <- dim(pfpr)
extent <- 579070:(579070 + 50)
in_bounds <- function(bounds, extent) {
    x <- numeric(length(extent))
    y <- numeric(length(extent))
    total <- 0
    for (i in 1:length(extent)) {
        z <- extent[i]
        xy <- decode_hilbert(z)
        if (xy[1] <= bounds[1] & xy[2] <= bounds[2]) {
            x[i] <- xy[1]
            y[i] <- xy[2]
            total <- total + 1
        }
    }
    if (total > 0) {
        list(row = x[1:total], col = y[1:total])
    } else {
        list(row = NULL, col = NULL)
    }
}

smallest_block <- function(image, row, col) {
    rowmin <- min(row)
    rows <- max(row) - rowmin + 1
    colmin <- min(col)
    cols <- max(col) - colmin + 1
    data <- matrix(
        raster::getValuesBlock(image, row = rowmin, nrows = rows, col = colmin, ncols = cols),
        byrow = TRUE,
        nrow = rows
    )
    list(data = data, rowoffset = rowmin - 1, coloffset = colmin - 1)
}
ib <- in_bounds(bounds, extent)
sb <- smallest_block(pfpr, ib$row, ib$col)

zmax <- 4057942
non_na_data(pfpr, am, in_bounds(dim(pfpr), as.integer(zmax / 4):(as.integer(zmax / 4) + 50)))

find_non_zero <- function(pfpr) {
    xmax <- raster::nrow(pfpr)
    ymax <- raster::ncol(pfpr)
    value <- NA
    while (is.na(value)) {
        x <- sample.int(xmax, 1)
        y <- sample.int(ymax, 1)
        value <- raster::extract(pfpr, matrix(c(x, y), ncol = 2))
    }
    list(x = x, y = y, value = value)
}


ld <- load_all("pr2ar_mesh.csv", c(1, 1), 32)
acc <- ld$pfpr
a0 <- acc[["2000"]]
a1 <- acc[["2001"]]
a0[4, 8] <- 7.2
a1[13, 17] <- 24
aa <- array(c(a0, a1), dim = c(32, 32, 2))
aa[4, 8, 1] == a0[4, 8]
aa[13, 17, 2] == a1[13, 17]
together <- array(do.call(c, acc), dim = c(32, 32, length(acc)))

output <- list()
for (i in 1:3) {
    output[[i]] <- list(block = c(row = i, col = 2 * i))
}
vapply(
    output,
    function(one_result) {
        one_result$block
    },
    FUN.VALUE = numeric(2)
    )


args <- check_args(arg_parser(c(
    "--config=gen_scaled_ar/rc_kappa.toml",
    "--country=gmb",
    "--years=2010:2011",
    "--overwrite",
    "--cores=1"
    )))

library(rampdata)
library(raster)
base <- ramp_path("/globalrc/outputs/basicr/201109")
year <- 2018
for (name in c("alpha", "kappa", "eir", "vc", "rc")) {
    fn <- as.path(add_path(base, file = sprintf("%s_%d.tif", name, year)))
    var <- raster::raster(fn)
    grDevices::png(file = sprintf("%s_%d.png", name, year))
    plot(var, main = sprintf("%s %d", name, year), useRaster = TRUE)
    grDevices::dev.off()
}

input_list <- list(
    pfpr = array(numeric(3*2*3), dim = c(3,2,3)),
    am = array(numeric(3*2*3), dim = c(3,2,3))
)
input_list$pfpr[1,1,1] <- NA
input_list$am[2,1,3] <- NA
not_available <- lapply(input_list, function(check) is.na(check))
not_available
dim_keep <- length(dim(input_list[[1]]))
keep <- !rowSums(array(do.call(c, not_available), dim = c(3,2,3,2)), dims = dim_keep)
dim(keep)

ff <- array(numeric(3*2*4), dim = c(3,2,3))
ff[keep] <- 1
ff

ff <- ll ~ rnorm(N, mean = 3, sd = 0.5)
ff
is.expression(parse(text = "rnorm(N, mean = 3, sd = 0.5)"))
eval(parse(text = "rnorm(N, mean = 3, sd = 0.5)"), envir = list(N=5))
if (length(grep("(", "hi(", fixed = TRUE)>0)) cat(3)
grep("(", "hi(there", fixed = TRUE)

ll <- data.frame(a = c(1,2), b = c(3,4))
one_draw <- ll[1,]
one_draw
with(one_draw, cat(sprintf("a %f b %f\n", a, b)))

arr <- array(runif(8), dim = c(2,2,2))
dim(arr)
typeof(arr)
class(arr)
is.array(arr)
stats::quantile(arr, c(0.025, 0.5, .975), dim = c(1,2))
apply(arr, MARGIN = c(2, 3), FUN = function(x) stats::quantile(x, c(0.025, 0.5, 0.975)))
dim(arr) <- c(2, 4)

qgamma(0.5, c(4, 4), c(1, 2))

a <- if (3 > 2) {
    7
} else {
    9
}
ll <- list(
    a = list(lower = c(1,2), median = c(3,4), upper = c(5,6)),
    b = list(lower = c(1,2), median = c(3,4), upper = c(5,6)),
    c = list(lower = c(1,2), median = c(3,4), upper = c(5,6))
)
lapply(1:(length(ll[[1]])*length(ll)), function(idx) {
    list_idx <- 1 + (idx - 1) %/% length(ll)
    sublist_idx <- 1 + idx %% length(ll)
    name <- names(ll)[list_idx]
    sub_list <- ll[[list_idx]]
    sub_name <- names(sub_list)[sublist_idx]
    paste(name, sub_name, sep = "_")
})
ll[[1]]
seq(0, 9, 1)
seq(0, 9, 1) %/% 3
seq(0, 9, 1) %% 3


