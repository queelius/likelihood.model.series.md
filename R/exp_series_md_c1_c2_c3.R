#' Constructs a likelihood model for `exp_series_md_c1_c2_c3`.
#'
#' Likelihood model for exponential series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3,
#' described below.
#'
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood.model` package by providing the following methods:
#'
#'  (1) `loglik.exp_series_md_c1_c2_c3`
#'  (2) `score.exp_series_md_c1_c2_c3`
#'  (3) `hess_loglik.exp_series_md_c1_c2_c3`
#'
#' These are useful for doing maximum likelihood estimation, hypothesis
#' testing (e.g., likelihood ratio test), estimation of asymptotic sampling
#' distribution given data from the DGP according to the specified model,
#' etc.
#'
#' It is designed to work well with the `likelihood_model` R package. In
#' particular, it is intended to be used with the `likelihood_contr_model`
#' object, which is a `likelihood_model` object that allows likelihood
#' contributions to be added for whatever data model you have in mind.
#'
#' In this likelihood model, masked component data approximately satisfies the
#' following conditions:
#'
#' C1: `Pr{K[i] in C[i]) = 1`
#' C2: `Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]) = Pr(C[i]=c[i] | K[i]=j', T[i]=t[i])`
#'     for any `j, j' in c[i]`.
#' C3: masking probabilities are independent of `theta`
#'
#' As a special case, this model also includes exact component cause of failure
#' data where the candidate set is a singleton.
#'
#' @param rates rate parameters for exponential component lifetimes (optional,
#'   used as initial values for MLE if provided)
#' @param lifetime column name for system lifetime, defaults to `"t"`
#' @param indicator column name for right-censoring indicator, defaults to
#'   `"delta"`. TRUE/1 = exact failure time, FALSE/0 = right-censored.
#'   For backwards compatibility, if this column is not present in the data,
#'   censoring is inferred from empty candidate sets (all FALSE).
#' @param candset column prefix for candidate set indicators, defaults to `"x"`
#' @export
#' @return likelihood model object with class
#'   `c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")`
#' @examples
#' # Create model and fit to data using generic dispatch
#' model <- exp_series_md_c1_c2_c3()
#' # solver <- fit(model)
#' # mle <- solver(data, par = c(1, 1, 1))
exp_series_md_c1_c2_c3 <- function(rates = NULL, lifetime = "t",
                                   indicator = "delta", candset = "x") {
  structure(
    list(
      rates = rates,
      lifetime = lifetime,
      indicator = indicator,
      candset = candset
    ),
    class = c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")
  )
}


#' Log-likelihood method for `exp_series_md_c1_c2_c3` model.
#'
#' Returns a log-likelihood function for an exponential series system with
#' respect to rate parameters for masked data with candidate sets that satisfy
#' conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: rate parameters of exponential component lifetime distributions
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `indicator`: right-censoring indicator column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model loglik
#' @method loglik exp_series_md_c1_c2_c3
#' @export
loglik.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime, indicator = defaults$indicator,
           candset = defaults$candset, ...) {
    if (any(par <= 0)) return(-Inf)
    d <- extract_model_data(df, lifetime, indicator, candset)
    if (length(par) != d$m) {
      stop(sprintf("Expected %d parameters but got %d", d$m, length(par)))
    }

    # Log-likelihood: -sum(t)*sum(theta) + sum_{i: delta_i=TRUE} log(sum(theta[C_i]))
    ll <- -sum(d$t) * sum(par)
    for (i in seq_len(d$n)) {
      if (d$delta[i]) {
        if (!any(d$C[i, ])) {
          stop("C1 violated: exact observation with empty candidate set at row ", i)
        }
        theta_c <- sum(par[d$C[i, ]])
        if (theta_c > 0) {
          ll <- ll + log(theta_c)
        }
      }
    }
    ll
  }
}


#' Score method for `exp_series_md_c1_c2_c3` model.
#'
#' Returns a score (gradient) function for an exponential series system with
#' respect to parameter `theta` for masked component failure with candidate
#' sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return score function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: rate parameters of exponential component lifetime distributions
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `indicator`: right-censoring indicator column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model score
#' @method score exp_series_md_c1_c2_c3
#' @export
score.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime, indicator = defaults$indicator,
           candset = defaults$candset, ...) {
    if (any(par <= 0)) return(rep(NA, length(par)))
    d <- extract_model_data(df, lifetime, indicator, candset)
    if (length(par) != d$m) {
      stop(sprintf("Expected %d parameters but got %d", d$m, length(par)))
    }

    # Score: d/d(theta_j) = -sum(t) + sum_{i: delta_i} C[i,j] / (C[i,] %*% theta)
    v <- rep(-sum(d$t), d$m)
    if (any(d$delta)) {
      c_delta <- d$C[d$delta, , drop = FALSE]
      theta_c <- as.numeric(c_delta %*% par)
      v <- v + as.numeric(colSums(c_delta / theta_c))
    }
    v
  }
}


#' Hessian of log-likelihood method for `exp_series_md_c1_c2_c3` model.
#'
#' Returns the observed information matrix (negative Hessian) for an
#' exponential series system with respect to parameter `theta` for masked data
#' with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return hessian function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: rate parameters of exponential component lifetime distributions
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `indicator`: right-censoring indicator column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model hess_loglik
#' @method hess_loglik exp_series_md_c1_c2_c3
#' @export
hess_loglik.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime, indicator = defaults$indicator,
           candset = defaults$candset, ...) {
    if (any(par <= 0)) return(matrix(NA, length(par), length(par)))
    d <- extract_model_data(df, lifetime, indicator, candset)
    if (length(par) != d$m) {
      stop(sprintf("Expected %d parameters but got %d", d$m, length(par)))
    }

    # Hessian: H[j,k] = -sum_{i: delta_i} C[i,j]*C[i,k] / (C[i,] %*% theta)^2
    hess <- matrix(0, nrow = d$m, ncol = d$m)
    if (any(d$delta)) {
      c_delta <- d$C[d$delta, , drop = FALSE]
      theta_c <- as.numeric(c_delta %*% par)
      weights <- -1 / theta_c^2
      hess <- unname(crossprod(c_delta * weights, c_delta))
    }
    hess
  }
}


#' Assumptions for `exp_series_md_c1_c2_c3` model.
#'
#' Returns the assumptions made by this likelihood model.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (ignored)
#' @return character vector of assumptions
#' @importFrom likelihood.model assumptions
#' @method assumptions exp_series_md_c1_c2_c3
#' @export
assumptions.exp_series_md_c1_c2_c3 <- function(model, ...) {
  c(
    SERIES_SYSTEM_ASSUMPTIONS[1],
    "exponential component lifetimes",
    SERIES_SYSTEM_ASSUMPTIONS[2],
    MASKING_CONDITIONS
  )
}


#' Random data generation for `exp_series_md_c1_c2_c3` model.
#'
#' Returns a function that generates random masked data from the exponential
#' series system DGP at a given parameter value.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return function that takes (theta, n, tau, p, ...) and returns a data frame
#'   with columns: t, delta, x1, x2, ..., xm
#' @importFrom likelihood.model rdata
#' @importFrom stats rexp runif
#' @method rdata exp_series_md_c1_c2_c3
#' @export
rdata.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(theta, n, tau = Inf, p = 0, ...) {
    m <- length(theta)
    if (any(theta <= 0)) stop("All rates must be positive")
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

    comp_lifetimes <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_lifetimes[, j] <- rexp(n, rate = theta[j])
    }

    generate_masked_series_data(
      comp_lifetimes, n, m, tau, p,
      defaults$lifetime, defaults$indicator, defaults$candset
    )
  }
}
