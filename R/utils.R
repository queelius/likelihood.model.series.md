# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Extract model column name defaults
#'
#' Helper function to extract default column names from a likelihood model
#' object. Used by all model methods to avoid repeating the same pattern.
#'
#' @param model likelihood model object with lifetime, indicator, candset fields
#' @return list with lifetime, indicator, candset defaults
#' @keywords internal
extract_model_defaults <- function(model) {
  list(
    lifetime = model$lifetime %||% "t",
    indicator = model$indicator %||% "delta",
    candset = model$candset %||% "x"
  )
}

#' Create numeric Hessian function from score function via Jacobian
#'
#' Factory function for creating Hessian functions using numerical
#' differentiation of the score function. Used by Weibull models which compute
#' score analytically but need numerical Hessian.
#'
#' @param score_fn score function returned by score.* method
#' @param model likelihood model object (for extracting defaults)
#' @return function(df, par, ...) that computes the Hessian matrix
#' @importFrom numDeriv jacobian
#' @keywords internal
make_numeric_hessian <- function(score_fn, model) {
  defaults <- extract_model_defaults(model)
  function(df, par, lifetime = defaults$lifetime,
           indicator = defaults$indicator,
           candset = defaults$candset, ...) {
    # Use Richardson extrapolation of order 6 for numerical stability
    numDeriv::jacobian(
      func = function(theta) score_fn(df, theta, lifetime, indicator, candset),
      x = par,
      method.args = list(r = 6)
    )
  }
}

# Masking condition assumption strings (shared across all models)
# nolint start: object_name_linter
#' @keywords internal
MASKING_CONDITIONS <- c(
  "C1: failed component is in candidate set with probability 1",
  "C2: uniform probability for candidate sets given component cause",
  "C3: masking probabilities independent of system parameters"
)

#' @keywords internal
SERIES_SYSTEM_ASSUMPTIONS <- c(
  "iid observations",
  "series system configuration"
)
# nolint end

#' Extract and validate model data from a masked data frame
#'
#' Shared validation logic for all likelihood model methods. Checks that the
#' data frame is non-empty, the lifetime column exists, decodes the candidate
#' set matrix, and extracts the censoring indicator with backwards compatibility.
#'
#' @param df masked data frame
#' @param lifetime column name for system lifetime
#' @param indicator column name for right-censoring indicator
#' @param candset column prefix for candidate set indicators
#' @return list with components: t (lifetimes), delta (censoring indicators),
#'   C (candidate set matrix), m (number of components), n (number of obs)
#' @importFrom md.tools md_decode_matrix
#' @keywords internal
extract_model_data <- function(df, lifetime, indicator, candset) {
  n <- nrow(df)
  if (n == 0) stop("df is empty")
  if (!lifetime %in% colnames(df)) {
    stop("lifetime variable not in colnames(df)")
  }

  cmat <- md_decode_matrix(df, candset)
  if (is.null(cmat) || ncol(cmat) == 0) {
    stop("no candidate set found for candset prefix '", candset, "'")
  }
  m <- ncol(cmat)

  # Get censoring indicator (backwards compat: infer from candidate sets)
  if (indicator %in% colnames(df)) {
    delta <- as.logical(df[[indicator]])
  } else {
    delta <- rowSums(cmat) > 0
  }

  list(t = df[[lifetime]], delta = delta, C = cmat, m = m, n = n)
}

#' Generate masked series system data
#'
#' Shared data generation logic for all rdata methods. Takes pre-generated
#' component lifetimes and applies system lifetime calculation, right-censoring,
#' and candidate set generation.
#'
#' @param comp_lifetimes n x m matrix of component lifetimes
#' @param n number of observations
#' @param m number of components
#' @param tau right-censoring time
#' @param p masking probability for non-failed components
#' @param default_lifetime column name for system lifetime
#' @param default_indicator column name for censoring indicator
#' @param default_candset column prefix for candidate sets
#' @return data frame with system lifetime, censoring indicator, and candidate
#'   sets
#' @importFrom stats runif
#' @keywords internal
generate_masked_series_data <- function(comp_lifetimes, n, m, tau, p,
                                        default_lifetime, default_indicator,
                                        default_candset) {
  sys_lifetime <- apply(comp_lifetimes, 1, min)
  failed_comp <- apply(comp_lifetimes, 1, which.min)

  delta <- sys_lifetime <= tau
  sys_lifetime <- pmin(sys_lifetime, tau)

  # Candidate sets satisfying C1, C2, C3
  candset <- matrix(FALSE, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    if (delta[i]) {
      candset[i, failed_comp[i]] <- TRUE
      if (p > 0 && m > 1) {
        others <- seq_len(m)[-failed_comp[i]]
        candset[i, others] <- runif(length(others)) < p
      }
    }
  }

  df <- data.frame(sys_lifetime, delta)
  names(df) <- c(default_lifetime, default_indicator)

  for (j in seq_len(m)) {
    df[[paste0(default_candset, j)]] <- candset[, j]
  }

  df
}

#' Cumulative hazard function for a component hazard function
#'
#' Creates a cumulative hazard function from a hazard function by integrating.
#'
#' @param haz hazard function
#' @return A function that computes the cumulative hazard at time t
#' @importFrom stats integrate
#' @export
cum_haz <- function(haz) {
  function(t, ...) {
    integrate(haz, lower = 0, upper = t, ...)$value
  }
}

#' Quantile function for a component with custom survival function
#'
#' Finds the time t such that S(t) = p using root finding.
#' The survival function S(t) is assumed to be monotonically decreasing
#' from S(0) = 1 to S(inf) = 0.
#'
#' @param p probability (quantile level), must be in (0, 1)
#' @param surv survival function S(t, theta, ...)
#' @param theta parameter vector passed to surv
#' @param t_lower lower bound for search interval
#' @param t_upper upper bound for search interval (sqrt to avoid overflow)
#' @param ... additional arguments passed to surv
#' @return time t such that S(t) = p
#' @export
#' @examples
#' # Exponential survival function
#' surv_exp <- function(t, theta) exp(-theta * t)
#'
#' # Median lifetime (50th percentile) for rate = 2
#' qcomp(0.5, surv = surv_exp, theta = 2.0)
qcomp <- function(p, surv, theta,
                  t_lower = .Machine$double.eps,
                  t_upper = .Machine$double.xmax^0.5, ...) {
  stats::uniroot(
    f = function(t) surv(t, theta, ...) - p,
    lower = t_lower,
    upper = t_upper,
    tol = .Machine$double.eps^0.5
  )$root
}

#' Random generation for a component with custom survival function
#'
#' Generates random samples using inverse transform sampling.
#'
#' @param n number of samples to generate
#' @param surv survival function S(t, theta)
#' @param theta parameter vector passed to surv
#' @return vector of n random samples
#' @export
#' @examples
#' # Exponential survival function
#' surv_exp <- function(t, theta) exp(-theta * t)
#'
#' # Generate 10 random samples with rate = 2
#' set.seed(123)
#' rcomp(10, surv = surv_exp, theta = 2.0)
rcomp <- function(n, surv, theta) {
  p <- runif(n)
  sapply(p, qcomp, surv = surv, theta = theta)
}
