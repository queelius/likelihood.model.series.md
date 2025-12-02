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

qcomp <- function(p, haz, surv, theta, t0 = 1, ...) {
    stats::optim(
        par = t0,
        fn = function(t) surv(t, theta, ...) - p,
        gr = function(t) haz(t, theta, ...) / surv(t, theta, ...),
        method = "Brent",
        lower = 0,
        upper = Inf)$par
}

rcomp <- function(n, haz, surv, theta) {
    p <- runif(n)
    sapply(p, qcomp, haz = haz, surv = surv, theta = theta)
}
