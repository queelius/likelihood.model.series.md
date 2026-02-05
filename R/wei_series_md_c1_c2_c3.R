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
#' @param indicator column name for right-censoring indicator, defaults to
#'   `"delta"`. TRUE/1 = exact failure time, FALSE/0 = right-censored.
#'   For backwards compatibility, if this column is not present in the data,
#'   censoring is inferred from empty candidate sets (all FALSE).
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
                                   indicator = "delta", candset = "x") {
  structure(
    list(
      shapes = shapes,
      scales = scales,
      lifetime = lifetime,
      indicator = indicator,
      candset = candset
    ),
    class = c("wei_series_md_c1_c2_c3", "series_md", "likelihood_model")
  )
}


#' Log-likelihood method for `wei_series_md_c1_c2_c3` model.
#'
#' Returns a log-likelihood function for a Weibull series system with
#' respect to parameter vector (shape_1, scale_1, ..., shape_m, scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `indicator`: right-censoring indicator column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model loglik
#' @method loglik wei_series_md_c1_c2_c3
#' @export
loglik.wei_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime, indicator = defaults$indicator,
           candset = defaults$candset, ...) {
    d <- extract_model_data(df, lifetime, indicator, candset)

    k <- length(par)
    if (k != 2 * d$m) {
      stop(sprintf("Expected %d parameters but got %d", 2 * d$m, k))
    }

    shapes <- par[seq(1, k, 2)]
    scales <- par[seq(2, k, 2)]

    if (any(shapes <= 0) || any(scales <= 0)) return(-Inf)

    ll <- 0
    for (i in seq_len(d$n)) {
      ll <- ll - sum((d$t[i] / scales)^shapes)

      if (d$delta[i]) {
        cind <- d$C[i, ]
        if (!any(cind)) {
          stop("C1 violated: exact observation with empty candidate set at row ", i)
        }
        hazard_sum <- sum(
          shapes[cind] / scales[cind] * (d$t[i] / scales[cind])^(shapes[cind] - 1)
        )
        if (hazard_sum > 0) {
          ll <- ll + log(hazard_sum)
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
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return score function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `indicator`: right-censoring indicator column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model score
#' @method score wei_series_md_c1_c2_c3
#' @export
score.wei_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime, indicator = defaults$indicator,
           candset = defaults$candset, ...) {
    d <- extract_model_data(df, lifetime, indicator, candset)

    k <- length(par)
    if (k != 2 * d$m) {
      stop(sprintf("Expected %d parameters but got %d", 2 * d$m, k))
    }

    shapes <- par[seq(1, k, 2)]
    scales <- par[seq(2, k, 2)]

    if (any(shapes <= 0) || any(scales <= 0)) return(rep(NA, k))

    shape_scores <- rep(0, d$m)
    scale_scores <- rep(0, d$m)

    for (i in seq_len(d$n)) {
      rt_term_shapes <- -(d$t[i] / scales)^shapes * log(d$t[i] / scales)
      rt_term_scales <- (shapes / scales) * (d$t[i] / scales)^shapes

      mask_term_shapes <- rep(0, d$m)
      mask_term_scales <- rep(0, d$m)

      if (d$delta[i]) {
        cind <- d$C[i, ]
        if (any(cind)) {
          denom <- sum(
            shapes[cind] / scales[cind] * (d$t[i] / scales[cind])^(shapes[cind] - 1)
          )

          if (denom > 0) {
            numer_shapes <- 1 / d$t[i] * (d$t[i] / scales[cind])^shapes[cind] *
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

    scr <- rep(0, k)
    scr[seq(1, k, 2)] <- shape_scores
    scr[seq(2, k, 2)] <- scale_scores
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
#'   - `indicator`: right-censoring indicator column name (default from model)
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
#' @return function that takes (theta, n, tau, p, ...) and returns a data frame
#'   with columns: t, delta, x1, x2, ..., xm
#' @importFrom likelihood.model rdata
#' @importFrom stats rweibull runif
#' @method rdata wei_series_md_c1_c2_c3
#' @export
rdata.wei_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(theta, n, tau = Inf, p = 0, ...) {
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
      defaults$lifetime, defaults$indicator, defaults$candset
    )
  }
}
