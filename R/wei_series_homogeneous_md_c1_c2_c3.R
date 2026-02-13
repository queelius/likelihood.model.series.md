#' Constructs a likelihood model for `wei_series_homogeneous_md_c1_c2_c3`.
#'
#' Likelihood model for Weibull series systems with homogeneous shape parameter
#' and masked component cause of failure with candidate sets that satisfy
#' conditions C1, C2, and C3.
#'
#' This is a reduced model where all components share a common shape parameter k,
#' while retaining individual scale parameters. The parameter vector is
#' (k, scale_1, ..., scale_m), giving m+1 parameters instead of 2m.
#'
#' A key property of this model is that the series system lifetime is itself
#' Weibull distributed with shape k and scale lambda_s = (sum(scale_j^{-k}))^{-1/k}.
#'
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood.model` package by providing the following methods:
#'
#'  (1) `loglik.wei_series_homogeneous_md_c1_c2_c3`
#'  (2) `score.wei_series_homogeneous_md_c1_c2_c3`
#'  (3) `hess_loglik.wei_series_homogeneous_md_c1_c2_c3`
#'
#' In this likelihood model, masked component data approximately satisfies:
#'
#' C1: `Pr{K[i] in C[i]} = 1`
#' C2: `Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]} = Pr{C[i]=c[i] | K[i]=j', T[i]=t[i]}`
#'     for any `j, j' in c[i]`.
#' C3: masking probabilities are independent of `theta`
#'
#' @param shape common shape parameter for all Weibull component lifetimes
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
#'   `c("wei_series_homogeneous_md_c1_c2_c3", "series_md", "likelihood_model")`
#' @seealso [wei_series_md_c1_c2_c3()] for the full model with heterogeneous shapes
#' @examples
#' # Create model and fit to data using generic dispatch
#' model <- wei_series_homogeneous_md_c1_c2_c3()
#' # solver <- fit(model)
#' # theta: (shape, scale1, scale2, ...)
#' # mle <- solver(data, par = c(1.2, 1000, 900, 850))
wei_series_homogeneous_md_c1_c2_c3 <- function(shape = NULL, scales = NULL,
                                                lifetime = "t",
                                                lifetime_upper = "t_upper",
                                                omega = "omega",
                                                candset = "x") {
  structure(
    list(
      shape = shape,
      scales = scales,
      lifetime = lifetime,
      lifetime_upper = lifetime_upper,
      omega = omega,
      candset = candset
    ),
    class = c("wei_series_homogeneous_md_c1_c2_c3", "series_md",
              "likelihood_model")
  )
}


#' Log-likelihood method for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns a log-likelihood function for a Weibull series system with
#' homogeneous shape parameter. The parameter vector is (k, scale_1, ..., scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' Supports four observation types. Left-censored and interval-censored have
#' closed-form likelihood contributions because all shapes are equal.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: parameter vector (shape, scale1, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model loglik
#' @method loglik wei_series_homogeneous_md_c1_c2_c3
#' @export
loglik.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime,
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    d <- extract_model_data(df, lifetime, omega, candset, lifetime_upper)

    npars <- length(par)
    if (npars != d$m + 1) {
      stop(sprintf("Expected %d parameters but got %d", d$m + 1, npars))
    }

    k <- par[1]
    scales <- par[-1]

    if (k <= 0 || any(scales <= 0)) return(-Inf)

    # System parameters for homogeneous Weibull
    sum_beta_neg_k <- sum(scales^(-k))
    beta_sys <- sum_beta_neg_k^(-1 / k)

    ll <- 0
    for (i in seq_len(d$n)) {
      if (d$omega[i] %in% c("exact", "right")) {
        # Survival: exp(-sum((t/beta_j)^k)) = exp(-(t/beta_sys)^k)
        ll <- ll - sum((d$t[i] / scales)^k)

        if (d$omega[i] == "exact") {
          cind <- d$C[i, ]
          hazard_sum <- sum(
            k / scales[cind] * (d$t[i] / scales[cind])^(k - 1)
          )
          if (hazard_sum > 0) {
            ll <- ll + log(hazard_sum)
          }
        }
      } else if (d$omega[i] == "left") {
        # log(w_c) + log(1 - exp(-(tau/beta_sys)^k))
        cind <- d$C[i, ]
        w_c <- sum(scales[cind]^(-k)) / sum_beta_neg_k
        tau_i <- d$t[i]
        ll <- ll + log(w_c) + log(-expm1(-(tau_i / beta_sys)^k))
      } else if (d$omega[i] == "interval") {
        # log(w_c) + log(exp(-(a/beta_sys)^k) - exp(-(b/beta_sys)^k))
        cind <- d$C[i, ]
        w_c <- sum(scales[cind]^(-k)) / sum_beta_neg_k
        a <- d$t[i]
        b <- d$t_upper[i]
        u_a <- (a / beta_sys)^k
        u_b <- (b / beta_sys)^k
        # Stable: -u_a + log(1 - exp(-(u_b - u_a)))
        ll <- ll + log(w_c) - u_a + log(-expm1(-(u_b - u_a)))
      }
    }
    ll
  }
}


#' Score method for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns a score (gradient) function for a Weibull series system with
#' homogeneous shape parameter. The parameter vector is (k, scale_1, ..., scale_m)
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
#'   - `par`: parameter vector (shape, scale1, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model score
#' @importFrom numDeriv grad
#' @method score wei_series_homogeneous_md_c1_c2_c3
#' @export
score.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime,
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    d <- extract_model_data(df, lifetime, omega, candset, lifetime_upper)

    npars <- length(par)
    if (npars != d$m + 1) {
      stop(sprintf("Expected %d parameters but got %d", d$m + 1, npars))
    }

    k <- par[1]
    scales <- par[-1]

    if (k <= 0 || any(scales <= 0)) return(rep(NA, npars))

    # Analytical score for exact + right (existing formulas)
    shape_score <- 0
    scale_scores <- rep(0, d$m)

    for (i in seq_len(d$n)) {
      if (d$omega[i] %in% c("exact", "right")) {
        u <- d$t[i] / scales
        u_k <- u^k

        # d/dk [-sum(u_j^k)] = -sum(u_j^k * log(u_j))
        shape_score <- shape_score - sum(u_k * log(u))

        # d/d(scale_j) [-sum(u_j^k)] = k * u_j^k / scale_j
        scale_scores <- scale_scores + k * u_k / scales

        if (d$omega[i] == "exact") {
          cind <- d$C[i, ]
          if (any(cind)) {
            u_km1 <- u[cind]^(k - 1)
            h_j <- k / scales[cind] * u_km1
            sum_h <- sum(h_j)

            if (sum_h > 0) {
              # d(h_j)/dk = (1/scale_j) * u_j^(k-1) * (1 + k * log(u_j))
              dh_dk <- (1 / scales[cind]) * u_km1 * (1 + k * log(u[cind]))
              shape_score <- shape_score + sum(dh_dk) / sum_h

              # d(h_j)/d(scale_j) = -k * h_j / scale_j
              dh_dscale <- -k * h_j / scales[cind]
              scale_scores[cind] <- scale_scores[cind] + dh_dscale / sum_h
            }
          }
        }
      }
    }

    scr <- c(shape_score, scale_scores)

    # Numerical gradient for left + interval (closed-form loglik)
    li_indices <- which(d$omega %in% c("left", "interval"))
    if (length(li_indices) > 0) {
      ll_li <- function(theta) {
        kk <- theta[1]
        sc <- theta[-1]
        if (kk <= 0 || any(sc <= 0)) return(-Inf)

        sum_neg_k <- sum(sc^(-kk))
        beta_sys <- sum_neg_k^(-1 / kk)
        ll_val <- 0

        for (i in li_indices) {
          cind <- d$C[i, ]
          w_c <- sum(sc[cind]^(-kk)) / sum_neg_k

          if (d$omega[i] == "left") {
            tau_i <- d$t[i]
            ll_val <- ll_val + log(w_c) +
              log(-expm1(-(tau_i / beta_sys)^kk))
          } else {
            a <- d$t[i]
            b <- d$t_upper[i]
            u_a <- (a / beta_sys)^kk
            u_b <- (b / beta_sys)^kk
            ll_val <- ll_val + log(w_c) - u_a +
              log(-expm1(-(u_b - u_a)))
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


#' Hessian of log-likelihood method for `wei_series_homogeneous_md_c1_c2_c3`.
#'
#' Returns the Hessian (second derivative matrix) of the log-likelihood for a
#' Weibull series system with homogeneous shape. Computed numerically via the
#' Jacobian of the score.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return hessian function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: parameter vector (shape, scale1, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom numDeriv jacobian
#' @importFrom likelihood.model hess_loglik
#' @method hess_loglik wei_series_homogeneous_md_c1_c2_c3
#' @export
hess_loglik.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
  score_fn <- score.wei_series_homogeneous_md_c1_c2_c3(model, ...)
  make_numeric_hessian(score_fn, model)
}


#' Assumptions for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns the assumptions made by this likelihood model.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (ignored)
#' @return character vector of assumptions
#' @importFrom likelihood.model assumptions
#' @method assumptions wei_series_homogeneous_md_c1_c2_c3
#' @export
assumptions.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
  c(
    SERIES_SYSTEM_ASSUMPTIONS[1],
    "Weibull component lifetimes with COMMON shape parameter",
    SERIES_SYSTEM_ASSUMPTIONS[2],
    "system lifetime is Weibull with shape k and scale (sum(scale_j^{-k}))^{-1/k}",
    MASKING_CONDITIONS
  )
}


#' System scale parameter for homogeneous Weibull series
#'
#' For a series system with Weibull components sharing shape k but with
#' individual scales, the system lifetime is itself Weibull with shape k
#' and this computed scale.
#'
#' @param k common shape parameter
#' @param scales vector of component scale parameters
#' @return system scale parameter
#' @export
#' @examples
#' # 3-component system with common shape 1.2
#' wei_series_system_scale(k = 1.2, scales = c(1000, 900, 850))
wei_series_system_scale <- function(k, scales) {
  (sum(scales^(-k)))^(-1 / k)
}


#' Random data generation for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns a function that generates random masked data from the homogeneous
#' Weibull series system DGP at a given parameter value.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return function that takes (theta, n, tau, p, observe, ...) and returns a
#'   data frame with columns: t, omega, x1, x2, ..., xm
#' @importFrom likelihood.model rdata
#' @importFrom stats rweibull runif
#' @method rdata wei_series_homogeneous_md_c1_c2_c3
#' @export
rdata.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(theta, n, tau = Inf, p = 0, observe = NULL, ...) {
    if (length(theta) < 2) {
      stop("theta must have at least 2 elements (k, scale1, ...)")
    }
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

    k <- theta[1]
    scales <- theta[-1]
    m <- length(scales)

    if (k <= 0 || any(scales <= 0)) {
      stop("Shape and all scale parameters must be positive")
    }

    comp_lifetimes <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_lifetimes[, j] <- rweibull(n, shape = k, scale = scales[j])
    }

    generate_masked_series_data(
      comp_lifetimes, n, m, tau, p,
      defaults$lifetime, defaults$omega, defaults$candset,
      observe = observe
    )
  }
}


#' @method ncomponents wei_series_homogeneous_md_c1_c2_c3
#' @export
ncomponents.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
  if (is.null(model$scales)) NULL else length(model$scales)
}


#' @method component_hazard wei_series_homogeneous_md_c1_c2_c3
#' @export
component_hazard.wei_series_homogeneous_md_c1_c2_c3 <- function(model, j,
                                                                  ...) {
  function(t, par, ...) {
    k <- par[1]
    beta_j <- par[j + 1]
    (k / beta_j) * (t / beta_j)^(k - 1)
  }
}
