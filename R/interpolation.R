

#' Creates a function that interpolates an AR-to-PR dataset.
#' The data calculates attack rate as a function of PR and recovery.
#' @param pr_to_ar_dt A data.table with columns AR, PR, and rho.
#' @returns A function with signature function(pr, rho) -> AR.
#' @export
ar_of_pr_rho <- function(pr_to_ar_dt) {
  if (length(unique(pr_to_ar_dt$rho[1:4])) < 4) {
    # This is the route taken. flog.debug("PR moves faster in pr_to_ar mesh file.")
    dtrow <- sort(unique(pr_to_ar_dt$PR))
    dtcol <- sort(unique(pr_to_ar_dt$rho))
    stopifnot(min(dtrow) == 0)
    stopifnot(max(dtrow) == 1)
    stopifnot(min(dtcol) == 0)
    stopifnot(max(dtcol) == 1)
    stopifnot(nrow(pr_to_ar_dt) == length(dtrow) * length(dtcol))
    dtz <- array(pr_to_ar_dt$AR, dim = c(length(dtrow), length(dtcol)))
    function(pr, rho) {
      stopifnot(all(pr >= 0))
      stopifnot(all(pr <= 1))
      stopifnot(all(rho >= 0))
      stopifnot(all(rho <= 1))
      stopifnot(all(is.finite(pr)))
      stopifnot(all(is.finite(rho)))
      stopifnot(length(rho) > 0)  # This was the failure.
      stopifnot(length(pr) > 0)
      stopifnot(length(pr) == length(rho))
      # This call looks like x and y are mixed up, but check the tests below.
      tryCatch(
          pracma::interp2(dtcol, dtrow, dtz, xp = rho, yp = pr, method = "linear"),
	  error = function(e) {
              cat(paste("interp2 fail for rho", paste0(rho, collapse=",")))
              cat(paste("interp2 fail for pr", paste0(pr, collapse=",")))
	  }
      )
    }
  } else {
    flog.debug("rho moves faster in pr_to_ar mesh file.")
    dtrow <- sort(unique(pr_to_ar_dt$rho))
    dtcol <- sort(unique(pr_to_ar_dt$PR))
    stopifnot(min(dtrow) == 0)
    stopifnot(max(dtrow) == 1)
    stopifnot(min(dtcol) == 0)
    stopifnot(max(dtcol) == 1)
    stopifnot(nrow(pr_to_ar_dt) == length(dtrow) * length(dtcol))
    dtz <- array(pr_to_ar_dt$AR, dim = c(length(dtrow), length(dtcol)))
    function(pr, rho) {
      stopifnot(all(pr >= 0))
      stopifnot(all(pr <= 1))
      stopifnot(all(rho >= 0))
      stopifnot(all(rho <= 1))
      pracma::interp2(dtcol, dtrow, dtz, xp = pr, yp = rho, method = "linear")
    }
  }
}


# Keep this for testing the other one. This works for a single value.
ar_of_pr_rho2 <- function(pr_to_ar_dt) {
    dt <- pr_to_ar_dt
    function(pr, rho) {
        akima::interp(x = dt$rho, y = dt$PR, z = dt$AR, xo = rho, yo = pr,
        extrap = TRUE)[[3]]
    }
}


#' This goes the other way, from AR to PR, with rho=0.
#' @param pr_ar_data This is the table of pr, ar, and rho values.
#' @export
build_ar2pr <- function(pr_ar_data) {
  no_treatment <- pr_ar_data[rho < 1e-6, c("PR", "AR")]
  # The sample data doesn't go all the way to 0 and 1, but that's the asymptotic value.
  AR_sample <- c(0, no_treatment$AR, 1)
  PR_sample <- c(0, no_treatment$PR, 1)
  function(attack_rates) {
    pracma::interp1(AR_sample, PR_sample, attack_rates, method = "linear")
  }
}
