
kappa_rm <- function(pfpr, c) { c * pfpr }
alpha_from_eir <- function(eir, pfpr, rho) {
  alpha <- 1 - exp(-log(1 + eir * k * tau) / k)
  alpha + ar2pr(pfpr, rho) - ar2pr(pfpr, 0)
}
rc_basic <- function(b, c, r, V, am) { b * V * c / r }

#' Compute the results of analyzing a time series of PfPR and AM.
#' @param pfpr A numeric vector of PR. Length is number of time points.
#' @param am A numeric vector of case management measure. Length is number of time points.
#' @param pr2ar A function that gives AR from PR and rho.
#' @param params Other parameters for the calculation in a list.
#' @return A list where each variable has results for the time series
#'     as a vector of numeric. The entries are `alpha`, `kappa`, `eir`, `vc`, `rc`.
#'
#' The parameters passed into this function are defined in the config
#' file, whose default name is `rc_kappa.toml`, in the `[parameters]` section.
pixel_work <- function(pfpr, am, params, strategies) {
  # PfPR and AM come in with the same set of NA patterns, where there is no land.
  with(c(params, strategies), {
    alpha <- pr_to_ar(pfpr, kam * am)
    kappa <- kappaf(pfpr, c)
    h <- -log(1 - alpha) / tau
    eir <- (exp(h * k * tau) - 1) / (b * k * tau)
    V <- eir / kappa
    R <- rcf(b, c, r, V, am)
    list(
      alpha = alpha,
      kappa = kappa,
      eir = eir * 365,
      vc = V,
      rc = R
    )
  })
}


pr2eir=function(x, b=0.55, r=1/200, k=4.2){
  ((1-x)^-k-1)*(r/k/b)
}


#' Takes a parameter that is an expression (of drawing a distribution) and calls it.
draw_parameters <- function(parameters, N) {
    if (N > 1) {
        draw_params <- with(parameters, {
            data.frame(
                kam = rep(kam, N),
                r = rnorm(N, r, r_sd),
                k = k,
                ku = runif(N),  # Use as quantile within calculation.
                b = rbeta(N, b, b_shape1, b_shape2),
                c = rep(c, N),
                tau = rep(tau, N),
                D_low = rep(D_low, N),
                D_high = rep(D_high, N)
            )
        })
    } else {
        draw_params <- data.frame(
            kam = kam,
            r = r,
            k = k,  # mean for k calculation.
            ku = 0.5,  # Use as quantile within calculation.
            b = b,
            c = c,
            tau = tau,
            D_low = D_low,
            D_high = D_high
        )
    }
    # If any parameters are less than zero, we should redraw.
    # It's a rejection method to get truncated distributions.
    stopifnot(all(as.matrix(draw_params) > 0))
    draw_params
}


pixel_two <- function(pfpr, am, params, strategies) {
  # PfPR and AM come in with the same set of NA patterns, where there is no land.
  with(c(params, strategies), {
    rho <- kam * am
    # draw of q determined by previous uniform draw.
    kd <- qgamma(params$ku, k * (1 - pfpr)^.6 * 5, 5 * (1 - pfpr)^.6)
    stopifnot(length(kd) == length(pfpr))
    # We aren't using the mechanistic model to get an absolute value of alpha
    # because that alpha isn't close enough. We are using that model to
    # estimate how much treatment would shift PfPR, given an alpha.
    alpha_cronus <- pr_to_ar(pfpr, rho)
    # nt = no_treatment
    pfpr_nt <- ar2pr(alpha_cronus)
    # This EIR estimate is from a fit to data.
    max_deir <- 1500 / 365
    eir_nt <- pmin(
      pr2eir(pfpr_nt, b, r, kd),
      rep(Inf, length(pfpr_nt))
      )
    alpha_nt <- 1 - (eir_nt * b * kd * tau + 1)^(-1/kd)
    # Then go back and estimate how much treatment would mean alpha was higher.
    lower_ar <- pr_to_ar(pfpr, numeric(length(pfpr)))
    upper_ar <- pr_to_ar(pfpr, rho)
    # Move it the same relative distance to the upper bound.
    alpha <- 1 - (1 - alpha_nt) * (1 - upper_ar) / (1 - lower_ar)
    # alpha <- min(alpha, 1 - 1e-10)

    kappa <- c * pfpr_nt
    h <- -log(1 - alpha) / tau
    # alpha can be 1 after the shift, so cap FOI.
    # These will be infinite, and that's OK. We don't
    # get rid of NAN. Those are failures.
    max_finite <- 1e10
    h[h > max_finite] <- max_finite
    eir <- (exp(h * kd * tau) - 1) / (b * kd * tau)
    eir[eir > max_finite] <- max_finite
    V <- ifelse(
        kappa > 0,
        eir / kappa,
        0
    )
    V[V > max_finite] <- max_finite
    R <- b * V * c * ((1 - rho) * D_high + rho * D_low)
    R[R > max_finite] <- max_finite
    list(
      alpha = alpha,
      kappa = kappa,
      foi = h,
      # deir = daily eir, aeir = annual eir = 365 * deir.
      aeir = eir * 365,
      vc = V,
      rc = R
    )
  })
}
