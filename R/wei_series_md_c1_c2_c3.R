#' Constructs a likelihood model for `wei_series_md_c1_c2_c3`.
#'
#' Likelihood model for Weibull series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood.model` package by providing the following methods:
#'
#'  (1) `loglik.wei_series_md_c1_c2_c3`
#'  (2) `score.wei_series_md_c1_c2_c3`
#'  (3) `hess_loglik.wei_series_md_c1_c2_c3`
#'
#' The Weibull series system has 2m parameters:
#' (shape_1, scale_1, ..., shape_m, scale_m).
#'
#' In this likelihood model, masked component data approximately satisfies:
#'
#' C1: `Pr{K[i] in C[i]} = 1`
#' C2: `Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]} = Pr{C[i]=c[i] | K[i]=j', T[i]=t[i]}`
#'     for any `j, j' in c[i]`.
#' C3: masking probabilities are independent of `theta`
#'
#' @param shapes shape parameters for Weibull component lifetimes (optional)
#' @param scales scale parameters for Weibull component lifetimes (optional)
#' @param lifetime column name for system lifetime, defaults to `"t"`
#' @param lifetime_upper column name for interval upper bound, defaults to
#'   `"t_upper"`. Only used for interval-censored observations.
#' @param omega column name for observation type, defaults to `"omega"`.
#'   Must contain character values: `"exact"`, `"right"`, `"left"`, or
#'   `"interval"`.
#' @param candset column prefix for candidate set indicators, defaults to `"x"`
#' @export
#' @return likelihood model object with class
#'   `c("wei_series_md_c1_c2_c3", "series_md", "likelihood_model")`
#' @examples
#' # Create model and fit to data using generic dispatch
#' model <- wei_series_md_c1_c2_c3()
#' # solver <- fit(model)
#' # theta: (shape1, scale1, shape2, scale2, ...)
#' # mle <- solver(data, par = c(1, 1000, 1, 1000, 1, 1000))
wei_series_md_c1_c2_c3 <- function(shapes = NULL, scales = NULL, lifetime = "t",
                                   lifetime_upper = "t_upper",
                                   omega = "omega", candset = "x") {
  structure(
    list(
      shapes = shapes,
      scales = scales,
      lifetime = lifetime,
      lifetime_upper = lifetime_upper,
      omega = omega,
      candset = candset
    ),
    class = c("wei_series_md_c1_c2_c3", "series_md", "likelihood_model")
  )
}


#' Integrand for numerical integration of Weibull series likelihood
#'
#' Computes h_c(t) * S(t) where h_c = sum_{j in c} h_j(t) and
#' S(t) = exp(-sum_l H_l(t)). Used for left-censored and interval-censored
#' observations in the heterogeneous Weibull model.
#'
#' @param t time values (vector, for use with stats::integrate)
#' @param shapes shape parameters for all components
#' @param scales scale parameters for all components
#' @param cind logical vector indicating which components are in candidate set
#' @return vector of integrand values
#' @keywords internal
wei_series_integrand <- function(t, shapes, scales, cind) {
  sapply(t, function(ti) {
    if (ti <= 0) return(0)
    hc <- sum(
      shapes[cind] / scales[cind] * (ti / scales[cind])^(shapes[cind] - 1)
    )
    cum_h <- sum((ti / scales)^shapes)
    hc * exp(-cum_h)
  })
}


#' Log-likelihood method for `wei_series_md_c1_c2_c3` model.
#'
#' Returns a log-likelihood function for a Weibull series system with
#' respect to parameter vector (shape_1, scale_1, ..., shape_m, scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' Supports four observation types. Left-censored and interval-censored
#' observations require numerical integration (via stats::integrate) because
#' heterogeneous shapes prevent a closed-form solution.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom likelihood.model loglik
#' @importFrom stats integrate
#' @method loglik wei_series_md_c1_c2_c3
#' @export
loglik.wei_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime,
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    d <- extract_model_data(df, lifetime, omega, candset, lifetime_upper)

    k <- length(par)
    if (k != 2 * d$m) {
      stop(sprintf("Expected %d parameters but got %d", 2 * d$m, k))
    }

    shapes <- par[seq(1, k, 2)]
    scales <- par[seq(2, k, 2)]

    if (any(shapes <= 0) || any(scales <= 0)) return(-Inf)

    ll <- 0
    for (i in seq_len(d$n)) {
      if (d$omega[i] %in% c("exact", "right")) {
        ll <- ll - sum((d$t[i] / scales)^shapes)

        if (d$omega[i] == "exact") {
          cind <- d$C[i, ]
          hazard_sum <- sum(
            shapes[cind] / scales[cind] *
              (d$t[i] / scales[cind])^(shapes[cind] - 1)
          )
          if (hazard_sum > 0) {
            ll <- ll + log(hazard_sum)
          }
        }
      } else {
        # Left or interval: numerical integration of h_c(t)*S(t)
        cind <- d$C[i, ]

        if (d$omega[i] == "left") {
          lower <- if (any(shapes[cind] < 1)) .Machine$double.eps^0.5 else 0
          val <- stats::integrate(
            wei_series_integrand,
            lower = lower, upper = d$t[i],
            shapes = shapes, scales = scales, cind = cind,
            subdivisions = 200, rel.tol = 1e-8,
            stop.on.error = FALSE
          )$value
        } else {
          # interval
          val <- stats::integrate(
            wei_series_integrand,
            lower = d$t[i], upper = d$t_upper[i],
            shapes = shapes, scales = scales, cind = cind,
            subdivisions = 200, rel.tol = 1e-8,
            stop.on.error = FALSE
          )$value
        }

        if (val > 0) {
          ll <- ll + log(val)
        } else {
          return(-Inf)
        }
      }
    }
    ll
  }
}


#' Score method for `wei_series_md_c1_c2_c3` model.
#'
#' Returns a score (gradient) function for a Weibull series system with
#' respect to parameter vector (shape_1, scale_1, ..., shape_m, scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' Uses a hybrid approach: analytical score for exact and right-censored
#' observations, numerical gradient (via numDeriv) for left-censored and
#' interval-censored observations.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return score function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom likelihood.model score
#' @importFrom numDeriv grad
#' @method score wei_series_md_c1_c2_c3
#' @export
score.wei_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime,
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    d <- extract_model_data(df, lifetime, omega, candset, lifetime_upper)

    k <- length(par)
    if (k != 2 * d$m) {
      stop(sprintf("Expected %d parameters but got %d", 2 * d$m, k))
    }

    shapes <- par[seq(1, k, 2)]
    scales <- par[seq(2, k, 2)]

    if (any(shapes <= 0) || any(scales <= 0)) return(rep(NA, k))

    # Analytical score for exact + right
    shape_scores <- rep(0, d$m)
    scale_scores <- rep(0, d$m)

    for (i in seq_len(d$n)) {
      if (d$omega[i] %in% c("exact", "right")) {
        rt_term_shapes <- -(d$t[i] / scales)^shapes * log(d$t[i] / scales)
        rt_term_scales <- (shapes / scales) * (d$t[i] / scales)^shapes

        mask_term_shapes <- rep(0, d$m)
        mask_term_scales <- rep(0, d$m)

        if (d$omega[i] == "exact") {
          cind <- d$C[i, ]
          if (any(cind)) {
            denom <- sum(
              shapes[cind] / scales[cind] *
                (d$t[i] / scales[cind])^(shapes[cind] - 1)
            )

            if (denom > 0) {
              numer_shapes <- 1 / d$t[i] *
                (d$t[i] / scales[cind])^shapes[cind] *
                (1 + shapes[cind] * log(d$t[i] / scales[cind]))
              mask_term_shapes[cind] <- numer_shapes / denom

              numer_scales <- (shapes[cind] / scales[cind])^2 *
                (d$t[i] / scales[cind])^(shapes[cind] - 1)
              mask_term_scales[cind] <- numer_scales / denom
            }
          }
        }

        shape_scores <- shape_scores + rt_term_shapes + mask_term_shapes
        scale_scores <- scale_scores + rt_term_scales - mask_term_scales
      }
    }

    scr <- rep(0, k)
    scr[seq(1, k, 2)] <- shape_scores
    scr[seq(2, k, 2)] <- scale_scores

    # Numerical gradient for left + interval
    li_indices <- which(d$omega %in% c("left", "interval"))
    if (length(li_indices) > 0) {
      ll_li <- function(theta) {
        sh <- theta[seq(1, k, 2)]
        sc <- theta[seq(2, k, 2)]
        if (any(sh <= 0) || any(sc <= 0)) return(-Inf)

        ll_val <- 0
        for (i in li_indices) {
          cind <- d$C[i, ]

          if (d$omega[i] == "left") {
            lower <- if (any(sh[cind] < 1)) .Machine$double.eps^0.5 else 0
            val <- stats::integrate(
              wei_series_integrand,
              lower = lower, upper = d$t[i],
              shapes = sh, scales = sc, cind = cind,
              subdivisions = 200, rel.tol = 1e-8,
              stop.on.error = FALSE
            )$value
          } else {
            val <- stats::integrate(
              wei_series_integrand,
              lower = d$t[i], upper = d$t_upper[i],
              shapes = sh, scales = sc, cind = cind,
              subdivisions = 200, rel.tol = 1e-8,
              stop.on.error = FALSE
            )$value
          }

          if (val > 0) {
            ll_val <- ll_val + log(val)
          } else {
            return(-Inf)
          }
        }
        ll_val
      }

      grad_li <- numDeriv::grad(ll_li, par)
      scr <- scr + grad_li
    }

    scr
  }
}


#' Hessian of log-likelihood method for `wei_series_md_c1_c2_c3` model.
#'
#' Returns the Hessian (second derivative matrix) of the log-likelihood for a
#' Weibull series system. Computed numerically via the Jacobian of the score.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return hessian function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom numDeriv jacobian
#' @importFrom likelihood.model hess_loglik
#' @method hess_loglik wei_series_md_c1_c2_c3
#' @export
hess_loglik.wei_series_md_c1_c2_c3 <- function(model, ...) {
  score_fn <- score.wei_series_md_c1_c2_c3(model, ...)
  make_numeric_hessian(score_fn, model)
}


#' Assumptions for `wei_series_md_c1_c2_c3` model.
#'
#' Returns the assumptions made by this likelihood model.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (ignored)
#' @return character vector of assumptions
#' @importFrom likelihood.model assumptions
#' @method assumptions wei_series_md_c1_c2_c3
#' @export
assumptions.wei_series_md_c1_c2_c3 <- function(model, ...) {
  c(
    SERIES_SYSTEM_ASSUMPTIONS[1],
    "Weibull component lifetimes",
    SERIES_SYSTEM_ASSUMPTIONS[2],
    MASKING_CONDITIONS
  )
}


#' Random data generation for `wei_series_md_c1_c2_c3` model.
#'
#' Returns a function that generates random masked data from the Weibull
#' series system DGP at a given parameter value.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return function that takes (theta, n, tau, p, observe, ...) and returns a
#'   data frame with columns: t, omega, x1, x2, ..., xm
#' @importFrom likelihood.model rdata
#' @importFrom stats rweibull runif
#' @method rdata wei_series_md_c1_c2_c3
#' @export
rdata.wei_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(theta, n, tau = Inf, p = 0, observe = NULL, ...) {
    k <- length(theta)
    if (k %% 2 != 0) stop("theta must have even length (shape1, scale1, ...)")
    m <- k / 2
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

    shapes <- theta[seq(1, k, 2)]
    scales <- theta[seq(2, k, 2)]

    if (any(shapes <= 0) || any(scales <= 0)) {
      stop("All shape and scale parameters must be positive")
    }

    comp_lifetimes <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_lifetimes[, j] <- rweibull(n, shape = shapes[j], scale = scales[j])
    }

    generate_masked_series_data(
      comp_lifetimes, n, m, tau, p,
      defaults$lifetime, defaults$omega, defaults$candset,
      defaults$lifetime_upper, observe = observe
    )
  }
}


#' @method ncomponents wei_series_md_c1_c2_c3
#' @export
ncomponents.wei_series_md_c1_c2_c3 <- function(model, ...) {
  if (is.null(model$shapes)) NULL else length(model$shapes)
}


#' @method component_hazard wei_series_md_c1_c2_c3
#' @export
component_hazard.wei_series_md_c1_c2_c3 <- function(model, j, ...) {
  function(t, par, ...) {
    k_j <- par[2 * j - 1]
    beta_j <- par[2 * j]
    (k_j / beta_j) * (t / beta_j)^(k_j - 1)
  }
}
