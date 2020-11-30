
kappa_rm <- function(pfpr, c) { c * pfpr }

alpha_from_eir <- function(eir, pfpr, k, rho, tau) {
  alpha <- 1 - exp(-log(1 + eir * k * tau) / k)
  alpha + ar2pr(pfpr, rho) - ar2pr(pfpr, 0)
}
rc_basic <- function(b, c, r, V, am) { b * V * c / r }

pr2eir=function(x, b=0.55, r=1/200, k=4.2){
  ((1-x)^-k-1)*(r/k/b)
}

pr2k_bind = function(x, s1=6, s2 = 3, n=1){
  pr = rep(x, each=n)
  mn = 1.5 + s1*(1-pr)^1.3
  vr = .08 + s2*(1-pr)^1.3
  cbind(pr, pmax(1, stats::rnorm(length(pr), mn, vr)))
}


pr2eirS = function(x, n=1){
  pr = rep(x, each=n)
  r = stats::rnorm(length(pr), 1/200, 1/3000)
  k = pr2k_bind(pr)[, 2]
  b = stats::rbeta(length(pr), .55*100, .45*100)
  eir = pmin(((1-pr)^-k-1)*(r/k/b), 1500/365)
  cbind(pr=pr, deir=eir, aeir=eir*365, leir = log10(eir*365))
}


pr2k_quantile <- function(pfpr, ku, s1=6, s2 = 3) {
  mn <- 1.5 + s1 * (1 - pfpr)^1.3
  vr <- .08 + s2 * (1 - pfpr)^1.3
  pmax(1, stats::qnorm(ku, mn, vr))
}

pr2k1 = function(x, mx=1.8, s1=4, pw=1.6){
  mx + s1*(1-x)^pw
}

# There were two of these.
pr2k = function(x, fac=2.6){
  (pr2k1(x) + pr2k2(x))/fac
}

pr2k_three_quantile <- function(pfpr, ku, s1=6, s2 = 2) {
  mn <- pr2k(pfpr)
  vr <- .03 + s2 * (1 - pfpr)^1.3
  pmax(1, stats::qnorm(ku, mn, vr))
}

pr2deir_quantile <- function(pfpr, r, b, kd) {
  k <- pr2k_quantile(pfpr, kd)
  pmin(((1 - pfpr)^ - k - 1) * (r / k / b), 1500 / 365)
}

x2kappa = function(x, k1=.3, k2=0, pw=1){
  k1*x^pw *exp(-k2*x)
}


#' Takes a parameter that is an expression (of drawing a distribution) and calls it.
#' @param parameters The config file parameters, no draws on input.
#' @param N the number of draws to make. Can be 1 for no draws.
#' @return A data frame with one row per draw, even if it's one.
#' Goes with pixel_three.
#' @export
draw_parameters <- function(parameters, N) {
  if (N > 1) {
    draw_params <- with(parameters, {
      data.frame(
        kam = rep(kam, N),
        r_inv = stats::rnorm(N, 1 / r, r_inv_sd),
        k = k,
        ku = runif(N),  # Use as quantile within calculation.
        b = stats::rbeta(N, b * 100, (1 - b) * 100),
        c = stats::rbeta(N, c * 50, (1 - c) * 50),
        s = s,
        pfpr_min = rep(pfpr_min, N),
        pfpr_max = rep(pfpr_max, N)
      )
    })
  } else {
    draw_params <- data.frame(
      kam = kam,
      r_inv = 1 / r,
      k = k,  # mean for k calculation.
      ku = 0.5,  # Use as quantile within calculation.
      b = b,
      c = c,
      s = s,
      pfpr_min = pfpr_min,
      pfpr_max = pfpr_max
    )
  }
  # If any parameters are less than zero, we should redraw.
  # It's a rejection method to get truncated distributions.
  stopifnot(all(as.matrix(draw_params) > 0))
  draw_params
}


#' This version works harder on the low end and high end of PfPR
#' to make the values have good variance there.
#' @param pfpr A list of pfpr values.
#' @param am A list of treatment values, the same length as pfpr.
#' @param params either a list or data frame, to be used in a with-statement.
#' @param strategies A list of functions to call to do parts of the work.
#' @return a list with a member for each variable.
#' @export
pixel_four <- function(pfpr, am, params, strategies) {
  # PfPR and AM come in with the same set of NA patterns, where there is no land.
  with(c(params, strategies), {
    rho <- kam * am
    # We aren't using the mechanistic model to get an absolute value of alpha
    # because that alpha isn't close enough. We are using that model to
    # estimate how much treatment would shift PfPR, given an alpha.
    alpha_cronus <- pr_to_ar(pfpr, rho)
    # nt = no_treatment
    pfpr_nt <- ar2pr(alpha_cronus)
    pfpr_nt <- pmax(pmin(pfpr_nt, pfpr_max), pfpr_min)

    # draw of q determined by previous uniform draw.
    kd <- pr2k_three_quantile(pfpr_nt, ku)
    stopifnot(length(kd) == length(pfpr_nt))

    deir <- pr2eir(pfpr_nt, b, 1 / r_inv, kd)
    kappa <- x2kappa(pfpr_nt, k1 = c)
    V <- deir * (1 + s * kappa) / kappa

    D <- c * r_inv
    R <- b * V * D
    list(
      pfpr_nt = pfpr_nt,
      # deir = daily eir, aeir = annual eir = 365 * deir.
      aeir = deir * 365,
      alpha = alpha_cronus,
      kappa = kappa,
      rc = R
    )
  })
}


#' Takes a parameter that is an expression (of drawing a distribution) and calls it.
#' @param parameters The config file parameters, no draws on input.
#' @param N the number of draws to make. Can be 1 for no draws.
#' @return A data frame with one row per draw, even if it's one.
#' Goes with pixel_three.
#' @export
draw_parameters_three <- function(parameters, N) {
    if (N > 1) {
        draw_params <- with(parameters, {
            data.frame(
                kam = rep(kam, N),
                r = stats::rnorm(N, r, r_sd),
                k = k,
                ku = runif(N),  # Use as quantile within calculation.
                b = stats::rbeta(N, b, b_shape1, b_shape2),
                c = rep(c, N),
                tau = rep(tau, N),
                D_low = rep(D_low, N),
                D_high = rep(D_high, N),
                pfpr_min = rep(pfpr_min, N),
                pfpr_max = rep(pfpr_max, N)
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
            D_high = D_high,
            pfpr_min = pfpr_min,
            pfpr_max = pfpr_max
        )
    }
    # If any parameters are less than zero, we should redraw.
    # It's a rejection method to get truncated distributions.
    stopifnot(all(as.matrix(draw_params) > 0))
    draw_params
}


#' This version works harder on the low end and high end of PfPR
#' to make the values have good variance there.
#' @param pfpr A list of pfpr values.
#' @param am A list of treatment values, the same length as pfpr.
#' @param params either a list or data frame, to be used in a with-statement.
#' @param strategies A list of functions to call to do parts of the work.
#' @return a list with a member for each variable.
#' @export
pixel_three <- function(pfpr, am, params, strategies) {
  # PfPR and AM come in with the same set of NA patterns, where there is no land.
  with(c(params, strategies), {
    rho <- kam * am
    # draw of q determined by previous uniform draw.
    kd <- pr2k_quantile(pfpr, ku)
    stopifnot(length(kd) == length(pfpr))
    # We aren't using the mechanistic model to get an absolute value of alpha
    # because that alpha isn't close enough. We are using that model to
    # estimate how much treatment would shift PfPR, given an alpha.
    alpha_cronus <- pr_to_ar(pfpr, rho)
    # nt = no_treatment
    pfpr_nt <- ar2pr(alpha_cronus)
    pfpr_nt <- pmax(pmin(pfpr_nt, pfpr_max), pfpr_min)
    # This EIR estimate is from a fit to data.
    eir_nt <- pmin(((1 - pfpr_nt)^ - kd - 1) * (r / kd / b), 1500 / 365)
    alpha_nt <- 1 - (eir_nt * b * kd * tau + 1)^(-1/kd)
    # Then go back and estimate how much treatment would mean alpha was higher.
    lower_ar <- pr_to_ar(pfpr, numeric(length(pfpr)))
    upper_ar <- pr_to_ar(pfpr, rho)
    # If Cronus says to move AR 2/3 of the way to 1,
    # then move the AR we have 2/3 of the way to 1.
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
    R <- b * V * ((1 - rho) * D_high + rho * D_low)
    R[R > max_finite] <- max_finite
    list(
      arnotreat = lower_ar,
      artreat = upper_ar,
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


#' Takes a parameter that is an expression (of drawing a distribution) and calls it.
#' @param parameters A list of floats and integers.
#' @param N The integer number of parameters to draw. Can be 1 to use means.
draw_parameters_two <- function(parameters, N) {
    if (N > 1) {
        draw_params <- with(parameters, {
            data.frame(
                kam = rep(kam, N),
                r = stats::rnorm(N, r, r_sd),
                k = k,
                ku = runif(N),  # Use as quantile within calculation.
                b = stats::rbeta(N, b, b_shape1, b_shape2),
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

#' Central calculation for pixel work that has low Rc.
#' @param pfpr parasite rate.
#' @param am anti-malarial treatment rate.
#' @param params Parameters for the math.
#' @param strategies Functions to apply during calculation.
#' @return Vectors of results.
#' This was made with fits to the middle values of PfPR.
#' It was finished Monday 23 Nov 2020. This was the version sent
#' to MAP that had reasonable values (c is corrected in last step).
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
    # If Cronus says to move AR 2/3 of the way to 1,
    # then move the AR we have 2/3 of the way to 1.
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
    R <- b * V * ((1 - rho) * D_high + rho * D_low)
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
