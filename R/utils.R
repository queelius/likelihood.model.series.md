# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

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
#'         C (candidate set matrix), m (number of components), n (number of obs)
#' @importFrom md.tools md_decode_matrix
#' @keywords internal
extract_model_data <- function(df, lifetime, indicator, candset) {
    n <- nrow(df)
    if (n == 0) stop("df is empty")
    if (!lifetime %in% colnames(df)) stop("lifetime variable not in colnames(df)")

    C <- md_decode_matrix(df, candset)
    if (is.null(C) || ncol(C) == 0) {
        stop("no candidate set found for candset prefix '", candset, "'")
    }
    m <- ncol(C)

    # Get censoring indicator (backwards compatibility: if not present, infer from candidate sets)
    if (indicator %in% colnames(df)) {
        delta <- as.logical(df[[indicator]])
    } else {
        delta <- rowSums(C) > 0
    }

    list(t = df[[lifetime]], delta = delta, C = C, m = m, n = n)
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
#' @return data frame with system lifetime, censoring indicator, and candidate sets
#' @importFrom stats runif
#' @keywords internal
generate_masked_series_data <- function(comp_lifetimes, n, m, tau, p,
                                        default_lifetime, default_indicator,
                                        default_candset) {
    # System lifetime and failed component
    sys_lifetime <- apply(comp_lifetimes, 1, min)
    failed_comp <- apply(comp_lifetimes, 1, which.min)

    # Apply right-censoring
    delta <- sys_lifetime <= tau
    sys_lifetime <- pmin(sys_lifetime, tau)

    # Generate candidate sets satisfying C1, C2, C3
    # For exact observations: failed component always in set, others with prob p
    # For censored observations: empty candidate set
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

    # Build data frame
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

#' Quantile function for a component with custom hazard/survival
#'
#' Finds the time t such that S(t) = p using root finding.
#' The survival function S(t) is assumed to be monotonically decreasing
#' from S(0) = 1 to S(inf) = 0.
#'
#' @param p probability (quantile level), must be in (0, 1)
#' @param haz hazard function h(t, theta, ...) (currently unused but kept for API consistency)
#' @param surv survival function S(t, theta, ...)
#' @param theta parameter vector passed to surv
#' @param t_lower lower bound for search interval
#' @param t_upper upper bound for search interval (sqrt to avoid overflow in survival calculations)
#' @param ... additional arguments passed to surv
#' @return time t such that S(t) = p
#' @export
qcomp <- function(p, haz = NULL, surv, theta,
                  t_lower = .Machine$double.eps,
                  t_upper = .Machine$double.xmax^0.5, ...) {
    # Find t such that S(t) = p using root finding
    # S(t) is monotonically decreasing from S(0) = 1 to S(inf) = 0
    stats::uniroot(
        f = function(t) surv(t, theta, ...) - p,
        lower = t_lower,
        upper = t_upper,
        tol = .Machine$double.eps^0.5
    )$root
}

#' Random generation for a component with custom hazard/survival
#'
#' Generates random samples using inverse transform sampling.
#'
#' @param n number of samples to generate
#' @param haz hazard function h(t, theta)
#' @param surv survival function S(t, theta)
#' @param theta parameter vector passed to haz and surv
#' @return vector of n random samples
#' @export
rcomp <- function(n, haz, surv, theta) {
    p <- runif(n)
    sapply(p, qcomp, haz = haz, surv = surv, theta = theta)
}
