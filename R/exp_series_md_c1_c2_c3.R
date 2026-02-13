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
#' @param lifetime_upper column name for interval upper bound, defaults to
#'   `"t_upper"`. Only used for interval-censored observations.
#' @param omega column name for observation type, defaults to `"omega"`.
#'   Must contain character values: `"exact"` (failure at t), `"right"`
#'   (right-censored at t), `"left"` (left-censored: failed before t),
#'   or `"interval"` (failed in (t, t_upper)).
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
                                   lifetime_upper = "t_upper",
                                   omega = "omega", candset = "x") {
  structure(
    list(
      rates = rates,
      lifetime = lifetime,
      lifetime_upper = lifetime_upper,
      omega = omega,
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
#' Supports four observation types: exact failures, right-censored,
#' left-censored, and interval-censored. All have closed-form likelihood
#' contributions for the exponential model.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: rate parameters of exponential component lifetime distributions
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom likelihood.model loglik
#' @method loglik exp_series_md_c1_c2_c3
#' @export
loglik.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime,
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    if (any(par <= 0)) return(-Inf)
    d <- extract_model_data(df, lifetime, omega, candset, lifetime_upper)
    if (length(par) != d$m) {
      stop(sprintf("Expected %d parameters but got %d", d$m, length(par)))
    }

    lambda_sys <- sum(par)

    # Exact + right: survival contribution -lambda_sys * t_i
    er_idx <- d$omega %in% c("exact", "right")
    ll <- -sum(d$t[er_idx]) * lambda_sys

    # Exact: add hazard term log(lambda_c)
    for (i in which(d$omega == "exact")) {
      lambda_c <- sum(par[d$C[i, ]])
      ll <- ll + log(lambda_c)
    }

    # Left-censored: log(lambda_c) + log(1 - exp(-lambda_sys*tau)) - log(lambda_sys)
    for (i in which(d$omega == "left")) {
      lambda_c <- sum(par[d$C[i, ]])
      ll <- ll + log(lambda_c) + log(-expm1(-lambda_sys * d$t[i])) -
        log(lambda_sys)
    }

    # Interval-censored: log(lambda_c) - lambda_sys*a + log(1 - exp(-lambda_sys*(b-a))) - log(lambda_sys)
    for (i in which(d$omega == "interval")) {
      lambda_c <- sum(par[d$C[i, ]])
      a <- d$t[i]
      b <- d$t_upper[i]
      ll <- ll + log(lambda_c) - lambda_sys * a +
        log(-expm1(-lambda_sys * (b - a))) - log(lambda_sys)
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
#' All four observation types (exact, right, left, interval) have closed-form
#' score contributions.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return score function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: rate parameters of exponential component lifetime distributions
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom likelihood.model score
#' @method score exp_series_md_c1_c2_c3
#' @export
score.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime,
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    if (any(par <= 0)) return(rep(NA, length(par)))
    d <- extract_model_data(df, lifetime, omega, candset, lifetime_upper)
    if (length(par) != d$m) {
      stop(sprintf("Expected %d parameters but got %d", d$m, length(par)))
    }

    lambda_sys <- sum(par)

    # Exact + right: -sum(t_i) + sum_{i:exact} C[i,j] / lambda_c_i
    er_idx <- d$omega %in% c("exact", "right")
    v <- rep(-sum(d$t[er_idx]), d$m)

    exact_idx <- d$omega == "exact"
    if (any(exact_idx)) {
      c_exact <- d$C[exact_idx, , drop = FALSE]
      theta_c <- as.numeric(c_exact %*% par)
      v <- v + as.numeric(colSums(c_exact / theta_c))
    }

    # Left-censored: C[i,j]/lambda_c + tau*exp(-lambda_sys*tau)/(1-exp(-lambda_sys*tau)) - 1/lambda_sys
    for (i in which(d$omega == "left")) {
      cind <- d$C[i, ]
      lambda_c <- sum(par[cind])
      tau_i <- d$t[i]
      x <- lambda_sys * tau_i
      # tau*exp(-x)/(1-exp(-x)) = tau/(exp(x)-1)
      ratio <- tau_i / expm1(x)
      v <- v + as.numeric(cind) / lambda_c + ratio - 1 / lambda_sys
    }

    # Interval-censored: C[i,j]/lambda_c + N/D - 1/lambda_sys
    # where D = exp(-lambda*a) - exp(-lambda*b), N = -a*exp(-lambda*a) + b*exp(-lambda*b)
    for (i in which(d$omega == "interval")) {
      cind <- d$C[i, ]
      lambda_c <- sum(par[cind])
      a <- d$t[i]
      b <- d$t_upper[i]
      ea <- exp(-lambda_sys * a)
      eb <- exp(-lambda_sys * b)
      D <- ea - eb
      N <- -a * ea + b * eb
      v <- v + as.numeric(cind) / lambda_c + N / D - 1 / lambda_sys
    }

    v
  }
}


#' Hessian of log-likelihood method for `exp_series_md_c1_c2_c3` model.
#'
#' Returns the Hessian (second derivative matrix) of the log-likelihood for an
#' exponential series system with respect to parameter `theta` for masked data
#' with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' All four observation types have closed-form Hessian contributions.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return hessian function that takes the following arguments:
#'   - `df`: masked data frame
#'   - `par`: rate parameters of exponential component lifetime distributions
#'   - `lifetime`: system lifetime column name (default from model)
#'   - `lifetime_upper`: interval upper bound column name (default from model)
#'   - `omega`: observation type column name (default from model)
#'   - `candset`: prefix of Boolean matrix encoding candidate sets
#' @importFrom likelihood.model hess_loglik
#' @method hess_loglik exp_series_md_c1_c2_c3
#' @export
hess_loglik.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(df, par, lifetime = defaults$lifetime,
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    if (any(par <= 0)) return(matrix(NA, length(par), length(par)))
    d <- extract_model_data(df, lifetime, omega, candset, lifetime_upper)
    if (length(par) != d$m) {
      stop(sprintf("Expected %d parameters but got %d", d$m, length(par)))
    }

    lambda_sys <- sum(par)

    # Candidate-set contribution: -C[i,j]*C[i,k]/lambda_c^2
    # for exact, left, and interval observations
    has_cand_idx <- d$omega %in% c("exact", "left", "interval")
    hess <- matrix(0, nrow = d$m, ncol = d$m)
    if (any(has_cand_idx)) {
      c_cand <- d$C[has_cand_idx, , drop = FALSE]
      theta_c <- as.numeric(c_cand %*% par)
      weights <- -1 / theta_c^2
      hess <- unname(crossprod(c_cand * weights, c_cand))
    }

    # Scalar contributions (same value added to ALL entries of H)
    # Left: -tau^2*exp(-x)/(1-exp(-x))^2 + 1/lambda_sys^2
    # Interval: (N'*D - N^2)/D^2 + 1/lambda_sys^2
    scalar_sum <- 0

    for (i in which(d$omega == "left")) {
      tau_i <- d$t[i]
      x <- lambda_sys * tau_i
      ex <- exp(-x)
      denom <- 1 - ex
      scalar_sum <- scalar_sum - tau_i^2 * ex / denom^2 + 1 / lambda_sys^2
    }

    for (i in which(d$omega == "interval")) {
      a <- d$t[i]
      b <- d$t_upper[i]
      ea <- exp(-lambda_sys * a)
      eb <- exp(-lambda_sys * b)
      D <- ea - eb
      N <- -a * ea + b * eb
      Np <- a^2 * ea - b^2 * eb
      scalar_sum <- scalar_sum + (Np * D - N^2) / D^2 + 1 / lambda_sys^2
    }

    # Adding a scalar to a matrix adds it to every element in R
    hess <- hess + scalar_sum

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
#' @return function that takes (theta, n, tau, p, observe, ...) and returns a
#'   data frame with columns: t, omega, x1, x2, ..., xm
#' @importFrom likelihood.model rdata
#' @importFrom stats rexp runif
#' @method rdata exp_series_md_c1_c2_c3
#' @export
rdata.exp_series_md_c1_c2_c3 <- function(model, ...) {
  defaults <- extract_model_defaults(model)

  function(theta, n, tau = Inf, p = 0, observe = NULL, ...) {
    m <- length(theta)
    if (any(theta <= 0)) stop("All rates must be positive")
    if (p < 0 || p > 1) stop("p must be in [0, 1]")

    comp_lifetimes <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_lifetimes[, j] <- rexp(n, rate = theta[j])
    }

    generate_masked_series_data(
      comp_lifetimes, n, m, tau, p,
      defaults$lifetime, defaults$omega, defaults$candset,
      defaults$lifetime_upper, observe = observe
    )
  }
}


#' @method ncomponents exp_series_md_c1_c2_c3
#' @export
ncomponents.exp_series_md_c1_c2_c3 <- function(model, ...) {
  if (is.null(model$rates)) NULL else length(model$rates)
}


#' @method component_hazard exp_series_md_c1_c2_c3
#' @export
component_hazard.exp_series_md_c1_c2_c3 <- function(model, j, ...) {
  function(t, par, ...) {
    rep(par[j], length(t))
  }
}


#' @method conditional_cause_probability exp_series_md_c1_c2_c3
#' @export
conditional_cause_probability.exp_series_md_c1_c2_c3 <- function(model, ...) {
  function(t, par, ...) {
    probs <- par / sum(par)
    matrix(probs, nrow = length(t), ncol = length(par), byrow = TRUE)
  }
}


#' @method cause_probability exp_series_md_c1_c2_c3
#' @export
cause_probability.exp_series_md_c1_c2_c3 <- function(model, ...) {
  function(par, ...) {
    par / sum(par)
  }
}
