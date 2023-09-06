#' Exponential series systems
#'
#' This file contains functions related to the Exponential series distribution.
#' Functions include simulation, pdf, cdf, quantile, and other related
#' functions for the Exponential series distribution.
#' 
#' @name Exponential series system
#' @author Alex Towell
#' @keywords exponential, distribution, series, statistics
NULL

#' Qunatile function for exponential series.
#'
#' @param p Vector of quantiles.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param lower.tail Logical, if TRUE (default), probabilities are `P[X<=x]`
#'                   otherwise, `P[X > x]`.
#' @param log.p Logical, if TRUE, vector of probabilities `p` are given as `log(p)`.
#' @return Quantiles corresponding to the given probabilities `p`.
#' @importFrom stats qexp
#' @export
qexp_series <- function(p, rates, lower.tail = TRUE, log.p = FALSE) {
    qexp(p, sum(rates), lower.tail, log.p)
}

#' Survival function for exponential series.
#'
#' @param n Integer, number of observations to generate.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @return A vector of random variates from the specified exponential series
#'         distribution.
#' @importFrom stats rexp
#' @export
rexp_series <- function(n, rates, keep_latent = FALSE) {

    m <- length(rates)

    data <- NULL
    if (keep_latent) {
        data <- matrix(NA, nrow = n, ncol = m)
        for (j in 1:m) {
            data[, j] <- rexp(n, rates[j])
        }
        data <- cbind(apply(data, 1, min), data)
    } else {
        data <- rexp(n, sum(rates))
    }
    data
}

#' pdf for exponential series.
#'
#' @param t series system lifetime
#' @param rates rate parameters for exponential component lifetimes
#' @param log return the log of the pdf
#' @importFrom stats dexp
#' @export
dexp_series <- function(t, rates, log = FALSE) {
    dexp(t, sum(rates), log)
}

#' Cumulative distribution function for exponential series.
#'
#' @param t Vector of series system lifetimes.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param lower.tail Logical; if TRUE (default), probabilities are `P[X <= x],
#'                   otherwise, `P[X > x]`.
#' @param log.p Logical; if TRUE, return the log of the cdf.
#' @return The cumulative probabilities evaluated at the specified lifetimes.
#' @importFrom stats pexp
#' @export
pexp_series <- function(t, rates, lower.tail = TRUE, log.p = FALSE) {
    pexp(t, sum(rates), lower.tail, log.p)
}

#' Hazard function for exponential series.
#'
#' @param t Vector of series system lifetimes.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param log.p Logical; if TRUE, return the log of the hazard function.
#' @return The hazard function evaluated at the specified lifetimes.
#' @export
hazard_exp_series <- function(t, rates, log.p = FALSE) {
    ifelse(log.p,
        ifelse(t < 0, -Inf, sum(rates)),
        ifelse(t < 0, 0, sum(rates))
    )
}


#' Survival function for exponential series.
#'
#' @param t Vector of series system lifetimes.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param log.p Logical; if TRUE, return the log of the survival function.
#' @return The survival function evaluated at the specified lifetimes.
#' @importFrom stats pexp
#' @export
surv.exp_series <- function(t, rates, log.p = FALSE) {
    pexp(t, sum(rates), TRUE, log.p)
}

#' Mean function for exponential series.
#' 
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @return The mean of the exponential series distribution.
#' @export
#' @importFrom stats mean
mean.exp_series <- function(rates) {
    mean(rates)
}
