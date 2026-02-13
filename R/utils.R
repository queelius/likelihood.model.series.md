# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Extract model column name defaults
#'
#' Helper function to extract default column names from a likelihood model
#' object. Used by all model methods to avoid repeating the same pattern.
#'
#' @param model likelihood model object with lifetime, lifetime_upper, omega,
#'   candset fields
#' @return list with lifetime, lifetime_upper, omega, candset defaults
#' @keywords internal
extract_model_defaults <- function(model) {
  list(
    lifetime = model$lifetime %||% "t",
    lifetime_upper = model$lifetime_upper %||% "t_upper",
    omega = model$omega %||% "omega",
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
           lifetime_upper = defaults$lifetime_upper,
           omega = defaults$omega,
           candset = defaults$candset, ...) {
    # Use Richardson extrapolation of order 6 for numerical stability
    numDeriv::jacobian(
      func = function(theta) {
        score_fn(df, theta, lifetime, lifetime_upper, omega, candset)
      },
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
#' data frame is non-empty, required columns exist, decodes the candidate set
#' matrix, and validates observation types.
#'
#' @param df masked data frame
#' @param lifetime column name for system lifetime
#' @param omega column name for observation type. Must contain character values:
#'   "exact", "right", "left", or "interval".
#' @param candset column prefix for candidate set indicators
#' @param lifetime_upper column name for interval upper bound (required when
#'   interval-censored observations are present)
#' @return list with components: t (lifetimes), omega (character vector of
#'   observation types), C (candidate set matrix), m (number of components),
#'   n (number of observations), t_upper (upper bounds or NULL)
#' @keywords internal
extract_model_data <- function(df, lifetime, omega, candset,
                               lifetime_upper = NULL) {
  n <- nrow(df)
  if (n == 0) stop("df is empty")
  if (!lifetime %in% colnames(df)) {
    stop("lifetime variable not in colnames(df)")
  }
  if (!omega %in% colnames(df)) {
    stop("omega variable '", omega, "' not in colnames(df)")
  }

  cmat <- md_decode_matrix(df, candset)
  if (is.null(cmat) || ncol(cmat) == 0) {
    stop("no candidate set found for candset prefix '", candset, "'")
  }
  m <- ncol(cmat)

  # Read and validate omega column
  omega_vals <- as.character(df[[omega]])
  valid_types <- c("exact", "right", "left", "interval")
  invalid <- setdiff(unique(omega_vals), valid_types)
  if (length(invalid) > 0) {
    stop(
      "invalid omega values: ", paste(invalid, collapse = ", "),
      ". Must be one of: ", paste(valid_types, collapse = ", ")
    )
  }

  # Extract t_upper for interval-censored observations
  t_upper <- NULL
  if (!is.null(lifetime_upper) && lifetime_upper %in% colnames(df)) {
    t_upper <- df[[lifetime_upper]]
  }

  # Validate observations: all non-right observations need non-empty candidate sets
  for (i in seq_len(n)) {
    if (omega_vals[i] == "right") next

    if (!any(cmat[i, ])) {
      msg <- if (omega_vals[i] == "exact") {
        paste0("C1 violated: exact observation with empty candidate set at row ", i)
      } else {
        paste0(omega_vals[i], "-censored observation must have non-empty ",
               "candidate set at row ", i)
      }
      stop(msg)
    }

    if (omega_vals[i] == "interval") {
      if (is.null(t_upper)) {
        stop(
          "interval-censored observations require a '",
          lifetime_upper %||% "t_upper", "' column"
        )
      }
      if (t_upper[i] <= df[[lifetime]][i]) {
        stop(
          "interval-censored observation requires t_upper > t at row ", i
        )
      }
    }
  }

  list(
    t = df[[lifetime]], omega = omega_vals, C = cmat, m = m, n = n,
    t_upper = t_upper
  )
}

#' Generate masked series system data
#'
#' Shared data generation logic for all rdata methods. Takes pre-generated
#' component lifetimes and applies an observation mechanism, then generates
#' candidate sets satisfying conditions C1, C2, C3.
#'
#' @param comp_lifetimes n x m matrix of component lifetimes
#' @param n number of observations
#' @param m number of components
#' @param tau right-censoring time (used when \code{observe} is NULL)
#' @param p masking probability for non-failed components
#' @param default_lifetime column name for system lifetime
#' @param default_omega column name for observation type
#' @param default_candset column prefix for candidate sets
#' @param default_lifetime_upper column name for interval upper bound
#' @param observe observation functor created by \code{observe_*} functions.
#'   When NULL, uses \code{\link{observe_right_censor}(tau)} for backwards
#'   compatibility.
#' @return data frame with system lifetime, observation type, and candidate
#'   sets
#' @importFrom stats runif
#' @keywords internal
generate_masked_series_data <- function(comp_lifetimes, n, m, tau, p,
                                        default_lifetime, default_omega,
                                        default_candset,
                                        default_lifetime_upper = paste0(default_lifetime, "_upper"),
                                        observe = NULL) {
  sys_lifetime <- apply(comp_lifetimes, 1, min)
  failed_comp <- apply(comp_lifetimes, 1, which.min)

  if (is.null(observe)) observe <- observe_right_censor(tau)

  # Apply observation mechanism
  obs_t <- numeric(n)
  omega_vals <- character(n)
  t_upper <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    obs <- observe(sys_lifetime[i])
    obs_t[i] <- obs$t
    omega_vals[i] <- obs$omega
    t_upper[i] <- obs$t_upper
  }

  # Candidate sets satisfying C1, C2, C3
  # Apply masking to all non-right observations (exact, left, interval)
  candset <- matrix(FALSE, nrow = n, ncol = m)
  has_failure <- omega_vals != "right"
  for (i in which(has_failure)) {
    candset[i, failed_comp[i]] <- TRUE  # C1
    if (p > 0 && m > 1) {
      others <- seq_len(m)[-failed_comp[i]]
      candset[i, others] <- runif(length(others)) < p
    }
  }

  # Build data frame
  df <- data.frame(obs_t, omega_vals, stringsAsFactors = FALSE)
  names(df) <- c(default_lifetime, default_omega)

  # Add t_upper column only if any interval observations exist
  if (any(omega_vals == "interval")) {
    df[[default_lifetime_upper]] <- t_upper
  }

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
