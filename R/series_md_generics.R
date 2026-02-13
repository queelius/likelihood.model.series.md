#' Number of components in a series system model
#'
#' Returns the number of components m in the series system. If the model was
#' constructed without parameter values, returns NULL.
#'
#' @param model a likelihood model object
#' @param ... additional arguments (currently ignored)
#' @return integer number of components, or NULL if not determinable
#' @export
ncomponents <- function(model, ...) UseMethod("ncomponents")


#' Component hazard function
#'
#' Returns a closure computing the hazard function h_j(t; theta) for the j-th
#' component. The returned function takes the full parameter vector `par` and
#' extracts the relevant component parameters internally.
#'
#' @param model a likelihood model object
#' @param j component index (integer from 1 to m)
#' @param ... additional arguments passed to the returned closure (e.g.,
#'   covariates for proportional hazards extensions)
#' @return a function with signature `function(t, par, ...)` computing h_j(t)
#' @export
component_hazard <- function(model, j, ...) UseMethod("component_hazard")


#' Conditional cause-of-failure probability
#'
#' Returns a closure computing P(K=j | T=t, theta) for all components j,
#' conditional on a specific failure time t. By Theorem 6 of the foundational
#' paper, this equals h_j(t; theta) / sum_l h_l(t; theta).
#'
#' @param model a likelihood model object
#' @param ... additional arguments passed to the returned closure
#' @return a function with signature `function(t, par, ...)` returning an
#'   n x m matrix where n = length(t) and column j gives P(K=j | T=t, theta)
#' @export
conditional_cause_probability <- function(model, ...) {
  UseMethod("conditional_cause_probability")
}


#' Marginal cause-of-failure probability
#'
#' Returns a closure computing P(K=j | theta) for all components j,
#' marginalized over the system failure time T. By Theorem 5 of the foundational
#' paper, this equals E_T\[P(K=j | T, theta)\].
#'
#' The default method uses Monte Carlo integration via [rdata()].
#'
#' @param model a likelihood model object
#' @param ... additional arguments passed to the returned closure
#' @return a function with signature `function(par, ...)` returning an m-vector
#'   where element j gives P(K=j | theta)
#' @export
cause_probability <- function(model, ...) UseMethod("cause_probability")


# =============================================================================
# Default methods for series_md class (Tier 2 — derived from Tier 1)
# =============================================================================

#' @describeIn conditional_cause_probability Default for series_md models via
#'   component hazard ratios (Theorem 6)
#' @method conditional_cause_probability series_md
#' @export
conditional_cause_probability.series_md <- function(model, ...) {
  m_model <- ncomponents(model)
  function(t, par, ...) {
    m <- m_model %||% stop(
      "ncomponents(model) is NULL; cannot determine m from model alone. ",
      "Provide a model constructed with parameter values."
    )
    n <- length(t)
    H <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      h_j <- component_hazard(model, j)
      H[, j] <- h_j(t, par, ...)
    }
    row_sums <- rowSums(H)
    H / row_sums
  }
}


#' @describeIn cause_probability Default for series_md models via Monte Carlo
#'   integration (Theorem 5)
#' @method cause_probability series_md
#' @export
cause_probability.series_md <- function(model, ...) {
  cond_fn <- conditional_cause_probability(model)
  rdata_fn <- rdata(model)

  function(par, n_mc = 10000, tau = Inf, p = 0, ...) {
    df <- rdata_fn(theta = par, n = n_mc, tau = tau, p = p)
    # Use only exact (uncensored) observations for the expectation
    defaults <- extract_model_defaults(model)
    omega_col <- defaults$omega
    if (omega_col %in% colnames(df)) {
      is_exact <- df[[omega_col]] == "exact"
    } else {
      is_exact <- rep(TRUE, nrow(df))
    }
    t_exact <- df[[defaults$lifetime]][is_exact]
    if (length(t_exact) == 0) {
      m <- ncomponents(model)
      return(rep(NA_real_, m))
    }
    probs <- cond_fn(t_exact, par, ...)
    colMeans(probs)
  }
}
