

.default.params <- list(
    kam = 0.6, b = 0.55, c = 0.17, k = 40, r = 0.005, tau = 10,
    ku = 0.5, D_low = 5, D_high = 40
    )
.default.strategies <- list()
# Make a fake mesh by pinning down the four corners.
.pr2ar.mesh <- data.table::data.table(
        PR = c(0.0, 0.0, 1.0, 1.0),
        rho = c(0.0, 1.0, 0.0, 1.0),
        AR = c(0.0, 0.1, 1.0, 1.0)
    )


# Keep this for testing the other one. This works for a single value.
ar_of_pr_rho2 <- function(pr_to_ar_dt) {
    dt <- pr_to_ar_dt
    function(pr, rho) {
        akima::interp(x = dt$rho, y = dt$PR, z = dt$AR, xo = rho, yo = pr,
        extrap = TRUE)[[3]]
    }
}

test_that("pixel_work works", {
    pfpr <- c(0.0, 0.1, 0.2, 0.3, 0.4)
    am <- c(0.66, 0.6, 0.4, 0.2, 0.0)
    pr2ar <- ar_of_pr_rho(.pr2ar.mesh)
    .default.strategies$pr_to_ar <- pr2ar
    stopifnot("rho" %in% names(.pr2ar.mesh))
    ar2pr_data <- data.table::data.table(
        rho = c(0, 0),
        AR = c(0.01, 0.9),
        PR = c(0.01, 0.9)
    )
    .default.strategies$ar2pr <- build_ar2pr(ar2pr_data)
    results <- pixel_two(pfpr, am, .default.params, .default.strategies)
    # There will be some variables.
    expect_gt(length(results), 2)
    for (var in names(results)) {
        expect_equal(length(results[[var]]), length(pfpr))
        expect_true(all(!is.na(results[[var]])))
    }
})


test_that("pr2eirS matches the deterministic one", {
    pr <- 0.15
    cnt <- 100
    pfpr <- rep(pr, cnt)
    deir_a <- pr2eirS(pfpr)[, 2]
    deir_b <- pr2deir_quantile(pfpr, runif(cnt), runif(cnt), runif(cnt))
    expect_equal(length(deir_a), length(deir_b))
})
