.default.toml <- system.file("testdata", "rc_kappa.toml", package = "globalrc")

.arg.tests <- list(
    list(a = c("--country=UGA"), b = list(country = "UGA")),
    list(a = c("--outvars=33"), b = list(outvars = "33")),
    list(a = c("--overwrite"), b = list(overwrite = TRUE)),
    list(a = c("--years=2000:2008"), b = list(years = 2000:2008))
)

withr::with_file("pr2ar_mesh.csv", {
    write.table(data.frame(), file = "pr2ar_mesh.csv")
    for (atidx in 1:length(.arg.tests)) {
        test_that(paste("arg handling works", atidx), {
            result <- check_args(arg_parser(c(
                sprintf("--config=%s", .default.toml),
                .arg.tests[[atidx]]$a)))
            expected <- .arg.tests[[atidx]]$b
            for (name in names(expected)) {
                expect_equal(result[[name]], expected[[name]])
            }
        })
    }
})


test_that("pixel_bounding_box works has whole boundary correct", {
    rr <- raster::raster(nrow = 100, ncol = 10, xmn = 1, xmx = 3, ymn = 2, ymx = 4)
    ans <- pixel_bounding_box(rr, c(xmin = 1, xmax = 3, ymin = 2, ymax = 4))
    expect_equal(ans, c(rmin = 1, rmax = 100, cmin = 1, cmax = 10))
})


#' If this runs in a situation without data, then skip these tests.
test_requires_data <- function() {
    testthat::skip_if_not(
        dir.exists(rampdata::as.path(rampdata::ramp_path("/globalrc/inputs/PfPR_medians"))))
}


test_that("available_data finds version", {
    test_requires_data()
    rampdata::initialize_workflow(.default.toml)
    data <- available_data("201029", NULL, NULL)
    expect_equal(data$years, 2000:2019)
    expect_equal(data$whole_input_extent["rmax"], c(rmax = 1741))
    expect_equal(data$whole_input_extent["cmax"], c(cmax = 1681))
    expect_equal(data$domain_extent["rmax"], c(rmax = 1741))
    expect_equal(data$domain_extent["cmax"], c(cmax = 1681))
})


test_that("available_data finds UGA outline", {
    data <- available_data("201029", "uga", NULL)
    expect_lt(unname(data$domain_extent["rmax"]), 1741)
    expect_lt(unname(data$domain_extent["cmax"]), 1681)
})


test_that("years_in_filenames picks up any filename",{
    fns <- c("dfb2000.tif", "haha32.tif", "blah2019_16.tif")
    by_year <- years_in_filenames(fns)
    expect_equal(names(by_year), as.character(c(2000, 2019)))
    expect_equal(unname(by_year), c(fns[1], fns[3]))
})

.default.options <- list(blocksize = 32L, single_tile_max = 2000L)

test_that("plan_domain_decomposition makes one tile for small", {
    plan <- plan_domain_decomposition(list(domain_extent = c(rmin = 1, rmax = 32, cmin = 1001, cmax = 1032)),
        .default.options
    )
    decomp <- plan$tiles
    expect_equal(decomp$tile_cnt, c(row = 1L, col = 1L))
})


test_that("plan_domain_decomposition makes enough tiles to cover", {
    plan <- plan_domain_decomposition(list(domain_extent = c(rmin = 1, rmax = 320, cmin = 1, cmax = 1032)),
        .default.options
    )
    decomp <- plan$tiles
    expect_equal(decomp$tile_cnt, c(row = 10L, col = 33L))
})


test_that("task_work finds all the tiles", {
    work <- task_work(c(row = 5L, col = 13L))
    expect_equal(dim(work)[1], 2L)
    expect_equal(dim(work)[2], 5L * 13L)
})


test_that("task_work finds single", {
    work <- task_work(c(row = 1L, col = 1L))
    expect_equal(dim(work)[1], 2L)
    expect_equal(dim(work)[2], 1L)
})


test_that("raster_extent_from_work calculates single correctly", {
    work <- task_work(c(row = 1L, col = 1L))
    blocksize <- c(row = 16, col = 32)
    ext <- raster_extent_from_work(work, blocksize)
    expect_equal(ext["rmin"], c(rmin = 1))
    expect_equal(ext["rmax"], c(rmax = 16))
    expect_equal(ext["cmin"], c(cmin = 1))
    expect_equal(ext["cmax"], c(cmax = 32))
})


test_that("raster_extent_from_work calculates larger correctly", {
    work <- task_work(c(row = 3L, col = 5L))
    blocksize <- c(row = 16, col = 32)
    ext <- raster_extent_from_work(work, blocksize)
    expect_equal(ext["rmin"], c(rmin = 1))
    expect_equal(ext["rmax"], c(rmax = 3*16))
    expect_equal(ext["cmin"], c(cmin = 1))
    expect_equal(ext["cmax"], c(cmax = 5*32))
})


test_that("load_data loads", {
    test_requires_data()
    defaults <- arg_parser(sprintf("--config=%s", .default.toml))
    domain_extent <- c(rmin = 1L, rmax = 32L, cmin = 65L, cmax = 128L)
    data <- load_data(defaults$config, defaults$pr2ar, domain_extent, 2000:2001)
    expect_equal(data$parameters$b, 0.55)
    expect_equal(length(data$pfpr), 2L)
    pfpr_dim <- dim(data$pfpr[["2000"]])
    expect_equal(pfpr_dim, c(32, 64))
    am_dim <- dim(data$am[["2001"]])
    expect_equal(am_dim, c(32, 64))
    expect_equal(data$offset, c(row = 1L, col = 65L))
})


make_pfpr <- function() {
    pfpr <- array(dim = c(3, 5*16, 4*32))
    for (ic in 1:(4*32)) {
        for (ir in 1:(5*16)) {
            pfpr[1, ir, ic] <- 2000 + encode_hilbert(ir, ic)
        }
    }
    pfpr[2, , ] <- pfpr[1, , ] + 1
    pfpr[3, , ] <- pfpr[1, , ] + 2
    pfpr
}


test_that("collect_and_permute handles these lists", {
    pfpr <- make_pfpr()
    pfpr_list <- list("2000" = pfpr[1,,], "2001" = pfpr[2,,], "2002" = pfpr[3,,])
    out <- collect_and_permute(pfpr_list)
    expect_equal(out, pfpr)
})


test_that("prepare_timeseries", {
    # Make data that we can match up later.
    # 3 years here, 2000:2002.
    pfpr <- make_pfpr()
    am <- -pfpr
    data <- list(
        parameters = list(a = .1, b = .2),
        pr_to_ar_dt = data.table::data.table(),
        offset = c(row = 1, col = 1),
        pfpr = list("2000" = pfpr[1,,], "2001" = pfpr[2,,], "2002" = pfpr[3,,]),
        am = list("2000" = am[1,,], "2001" = am[2,,], "2002" = am[3,,]),
        years <- 2000:2002
    )
    work <- matrix(c(1, 1, 1, 3, 5, 3), nrow = 2, dimnames = list(c("row", "col"), NULL))
    blocksize <- c(row = 16, col = 32)
    process_extent <- c(rmin = 1, rmax = 80, cmin = 1, cmax = 128)
    ts <- prepare_timeseries(data, work, process_extent, blocksize)
    expect_equal(length(ts$chunks), 3L)
    expect_equal(names(ts$chunks[[1]]), c("tile", "pfpr", "am"))
    # 3rd chunk is tile row=5, col=3. Second year is 2001.
    pfpr2001 <- ts$chunks[[3]]$pfpr[2, , ]
    expect_equal(2001 + encode_hilbert((5-1) * 16 + 1, (3-1) * 32 + 1), pfpr2001[1, 1])
    am2002 <- ts$chunks[[3]]$am[3, , ]
    expect_equal(2002 + encode_hilbert((5-1) * 16 + 1, (3-1) * 32 + 1), -am2002[1, 1])
})


.default.params <- list(
    kam = 0.6, b = 0.55, c = 0.17, k = 40, r = 0.005, tau = 10,
    ku = 0.5, D_low = 5, D_high = 40
    )
.default.strategies <- list(kappaf = kappa_rm, rcf = rc_basic)
# Make a fake mesh by pinning down the four corners.
.pr2ar.mesh <- data.table::data.table(
        PR = c(0.0, 0.0, 1.0, 1.0),
        rho = c(0.0, 1.0, 0.0, 1.0),
        AR = c(0.0, 0.1, 1.0, 1.0)
    )


test_that("combine_output puts the right results in the right places", {
    blocksize <- c(row = 2, col = 3)
    blocked <- list()
    for (tile_idx in 1:3) {
        blocked[[tile_idx]] <- list(
            ar = array(1:24, dim = c(4, blocksize)),
            rc = array(1:24, dim = c(4, blocksize)),
            array_names = c("ar", "rc"),
            block = c(row = tile_idx, col = 2 * tile_idx)
        )
    }
    result <- combine_output(blocked, blocksize)
    expect_equal(names(result), c("ar", "rc"))
})


test_that("chunks going to pixel calc have data", {
    test_requires_data()
    args <- check_args(arg_parser(c(
        paste0("--config=", .default.toml),
        "--country=uga",  # gmb = The Gambia, for a small one. uga = Uganda.
        "--years=2010:2011",
        "--overwrite",
        "--cores=4"
    )))
    # By this point, someone told us what to do, so load all the data.
    whole <- c(row = 1741L, col = 1681L)
    base <- setNames(as.integer(round(whole / 2)), names(whole))
    process_extent <- c(rmin = 1, rmax = 64, cmin = 1, cmax = 32)
    domain <- process_extent + unname(base[c("row", "row", "col", "col")])
    blocksize <- c(row = 16, col = 16)
    work <- task_work(c(row = 4, col = 2))

    flog.debug(paste("extent", paste(domain, collapse = ","), "names", paste(names(domain), collapse = ","), "\n"))
    data <- load_data(args$config, args$pr2ar, domain, args$years)
    pfpr <- data$pfpr[["2010"]]
    expect_equal(dim(pfpr), c(64, 32))
    expect_equal(sum(is.na(pfpr)), 13)
    timeseries_data <- prepare_timeseries(data, work, process_extent, blocksize)

    nonzero <- chunk_nonzero(timeseries_data$chunks)
    # The values shouldn't be the same.
    sums <- chunk_sums(timeseries_data$chunks)
    expect_gt(length(unique(sums[1,])), 1)
    expect_gt(length(unique(sums[2,])), 1)
})


test_that("ar_of_pr_rho fits the function it should", {
    test_requires_data()
    pr2ar_rp <- rampdata::workflow_path("pr2ar")
    pr_to_ar_dt <- data.table::fread(rampdata::as.path(pr2ar_rp))
    pr_to_ar <- ar_of_pr_rho(pr_to_ar_dt)
    # Putting a unit test right here. Run every time. Hah.
    for (test_line in c(5, 20, 42)) {
        egval <- pr_to_ar_dt[test_line, ]
        relerr <- (pr_to_ar(egval$PR, egval$rho) - egval$AR) / egval$AR
        expect_lt(relerr, 0.01)
    }
})


test_that("draw_parameters can be used in a with statement", {
    params <- list(r = 0.005, r_sd = 1/3000, b = 0.55, b_shape1 = 55,
    k = 4.2, kam = 0.6, pfpr_min = 0.02, pfpr_max = 0.98,
    b_shape2 = 45, c = 0.17, tau = 10, D_low = 5, D_high = 40)
    draws <- draw_parameters(params, 100)
    ans = with(draws[37, ], {
        r * k * b * tau
    })
    expect_gt(ans, 0)
})


test_that("summarize_draws produces correct number of outputs", {
    ci <- 95
    draw_cnt <- 12
    x <- rgamma(draw_cnt, 9, 2)
    tile <- c(10,15)
    cube <- c(2, tile)
    draws_list <- lapply(1:draw_cnt, function(draw_idx) {
        same_draws <- array(rep(x[draw_idx], prod(cube)), dim = cube)
        list(A = same_draws, B = same_draws)
    })
    results <- summarize_draws(draws_list, ci)
    expect_equal(names(results), c("A_lower", "A_median", "A_upper",
        "B_lower", "B_median", "B_upper"))
})


test_that("flatten_draws does flatten", {
    ll <- list(
        a = list(lower = c(1,2), median = c(3,4), upper = c(5,6)),
        b = list(lower = c(1,2), median = c(3,4), upper = c(5,6)),
        c = list(lower = c(1,2), median = c(3,4), upper = c(42, 37)),
        d = list(lower = c(1,2), median = c(3,4), upper = c(42, 3)),
        e = list(lower = c(1,2), median = c(3,4), upper = c(4, 37))
    )
    fout <- flatten_draws(ll)
    expect_true("b_median" %in% names(fout))
    expect_equal(length(fout), length(ll) * length(ll[[1]]))
    expect_equal(fout[["c_upper"]], c(42, 37))
})


test_that("flatten_draws does flatten for length shorter", {
    ll <- list(
        a = list(lower = c(1,2), median = c(3,4), upper = c(5,6)),
        e = list(lower = c(1,2), median = c(3,4), upper = c(4, 37))
    )
    fout <- flatten_draws(ll)
    expect_true("a_median" %in% names(fout))
    expect_equal(length(fout), length(ll) * length(ll[[1]]))
    expect_equal(fout[["e_upper"]], c(4, 37))
})

any_na <- function(x) {
 lapply(x, function(y) all(!is.na(y)))
}

test_that("linearized_work finds NA sees all true", {
    arr1 <- array(1:12, c(2,3,2))
    arr2 <- array(1:12, c(2,3,2))
    input_list <- list(a = arr1, b = arr2)
    res = linearized_work(input_list, any_na)
    expect_true(all(vapply(res, function(y) all(!is.na(y)), FUN.VALUE=logical(1))))
})


expect_length <- function(n) {
    function(input_list) {
        right_length <- vapply(input_list, function(x) {length(x) == n}, logical(1))
        stopifnot(all(right_length))
        input_list
    }
}


test_that("linearized_work finds NA sees all true", {
    arr1 <- array(1:12, c(2,3,2))
    arr2 <- array(1:12, c(2,3,2))
    arr1[1,3,2] <- NA
    input_list <- list(a = arr1, b = arr2)
    res = linearized_work(input_list, expect_length(11))
    expect_true(is.na(res[["b"]][1,3,2]))
})



test_that("linearized_work finds NA sees double NA", {
    arr1 <- array(1:12, c(2,3,2))
    arr2 <- array(1:12, c(2,3,2))
    arr1[1,3,2] <- NA
    arr2[1,3,2] <- NA
    arr2[1,1,2] <- NA
    input_list <- list(a = arr1, b = arr2)
    res = linearized_work(input_list, expect_length(10))
    expect_true(is.na(res[["b"]][1,3,2]))
})



test_that("linearized_work finds NA sees second NA", {
    arr1 <- array(1:12, c(2,3,2))
    arr2 <- array(1:12, c(2,3,2))
    arr2[1,3,2] <- NA
    input_list <- list(a = arr1, b = arr2)
    res = linearized_work(input_list, expect_length(11))
    expect_true(is.na(res[["a"]][1,3,2]))
    expect_true(is.na(res[["b"]][1,3,2]))
})



test_that("my_tasks gives array with dims", {
    work <- data.table::data.table(
        row = c(1,2,3,4),
        col = c(5,6,7,8)
    )
    res <- my_tasks(1, 2, work)
    expect_equal(res, matrix(c(1, 5, 2, 6), nrow = 2, dimnames = list(c("row", "col"))))
    res <- my_tasks(2, 2, work)
    expect_equal(res, matrix(c(3, 7, 4, 8), nrow = 2, dimnames = list(c("row", "col"))))
})



test_that("my_tasks sees extra tasks", {
    work <- data.table::data.table(
        row = 1:7,
        col = 8:14
    )
    res <- my_tasks(1, 4, work)
    expect_equal(res, matrix(c(1, 8, 2, 9), nrow = 2, dimnames = list(c("row", "col"))))
    res <- my_tasks(2, 4, work)
    expect_equal(res, matrix(c(3, 10, 4, 11), nrow = 2, dimnames = list(c("row", "col"))))
    res <- my_tasks(4, 4, work)
    expect_equal(res, matrix(c(7, 14), nrow = 2, dimnames = list(c("row", "col"))))
})


test_that("my_tasks does further uneven", {
    work <- data.table::data.table(
        row = 1:7,
        col = 8:14
    )
    res <- my_tasks(1, 3, work)
    expect_equal(res, matrix(c(1, 8, 2, 9, 3, 10), nrow = 2, dimnames = list(c("row", "col"))))
    res <- my_tasks(2, 3, work)
    expect_equal(res, matrix(c(4, 11, 5, 12), nrow = 2, dimnames = list(c("row", "col"))))
    res <- my_tasks(3, 3, work)
    expect_equal(res, matrix(c(6, 13, 7, 14), nrow = 2, dimnames = list(c("row", "col"))))
})


test_that("my_tasks ok with too few", {
    work <- data.table::data.table(
        row = 1:7,
        col = 8:14
    )
    res <- my_tasks(1, 9, work)
    expect_equal(res, matrix(c(1, 8), nrow = 2, dimnames = list(c("row", "col"))))
    res <- my_tasks(8, 9, work)
    expect_true(is.null(res))
    res <- my_tasks(9, 9, work)
    expect_true(is.null(res))
    res <- my_tasks(7, 9, work)
    expect_equal(res, matrix(c(7, 14), nrow = 2, dimnames = list(c("row", "col"))))
})


mt_cnt <- 5
row_cnt <- sample(1:100, mt_cnt, replace = TRUE)
task_cnt <- sample(1:100, mt_cnt, replace = TRUE)
for (mt_idx in 1:mt_cnt) {
    test_that(sprintf("my_tasks always gets all parts %d", mt_idx), {
        work <- data.table::data.table(
            row = 1:row_cnt[mt_idx],
            col = 1:row_cnt[mt_idx]
        )
        task_cnt <- task_cnt[mt_idx]
        task_set <- numeric(0)
        for (task_idx in 1:task_cnt) {
            ns <- my_tasks(task_idx, task_cnt, work)["row", ]
            expect_equal(length(intersect(task_set, ns)), 0)
            task_set <- c(task_set, ns)
        }
        expect_equal(sort(work$row), unname(sort(task_set)))
    })
}
