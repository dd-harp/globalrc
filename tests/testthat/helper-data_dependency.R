#' If this runs in a situation without data, then skip these tests.
test_requires_data <- function() {
    testthat::skip_if_not(
        tryCatch({
            dir.exists(
                rampdata::as.path(
                    rampdata::ramp_path("/globalrc/inputs/PfPR_medians")))
        },
        error = function(x) {
            FALSE
        }
    ))
}
