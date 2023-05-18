#' General series
#'
#' This file contains functions related to a general series distribution.
#' Functions include simulation, pdf, cdf, quantile, and other related
#' functions.
#' 
#' For parameter estimation from masked data, see the file \code{md_gen_series_mle.R}.
#'
#' @author Alex Towell
#' @name General series
#' @keywords distribution, series, statistics
NULL



#' Quantile function (inverse of the cdf) for a general series system.
#' By definition, the quantile \code{p * 100%} is the value \code{t} that
#' satisfies \code{F(t) - p = 0}. We solve for \code{t} using Newton's method.
#'
#' @param p vector of probabilities.
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of reliability functions, default is NULL (numerical approximation from h)
#' @param eps stopping condition, default is 1e-3
#' @param t0 initial guess, default is 1
#' @export
qgen_series <- Vectorize(function(p,theta,nparams,h,R=NULL,eps=1e-3,t0=1)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(eps > 0)
    stopifnot(all(p > 0))
    m <- length(h)

    if (is.null(R))
    {
        R <- list()
        for (i in 1:m)
            R[[i]] <- generate_survival_from_hazard(h[[i]])
    }
    stopifnot(length(h)==length(R))

    series_h <- hazard_gen_series_helper(h,nparams)
    series_R <- survival_gen_series_helper(R,nparams)
    t1 <- NULL
    repeat
    {
        alpha <- 1
        repeat
        {
            t1 <- t0 - alpha * (1-series_R(t,theta))/(series_h(t,theta)*series_R(t,theta))
            if (t1 > 0)
                break
            alpha <- alpha / 2
        }
        if (abs(t2-t1) < eps)
            break
        t1 <- t2
    }
    t2
}, vectorize.args = "p")

#' Sample from a general series system of \code{m} components whose hazard
#' and reliability functions are respectively given by \code{h} and \code{R}.
#'
#' @param n sample size
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of reliability functions, defaults to NULL (results in a numerical approximation from h)
#' @export
rgen_series <- function(n,theta,nparams,h,R=NULL)
{
    qgen_series(runif(n),theta,nparams,h,R)
}

#' pdf for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of reliability functions
#' @export
dgen_series <- Vectorize(function(t,theta,nparams,h,R)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    series_h <- hazard_gen_series_helper(h,nparams)
    series_R <- survival_gen_series_helper(R,nparams)
    series_h(t,theta) / series_R(t,theta)
}, vectorize.args="t")

#' cdf for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param R list of reliability functions
#' @export
pgen_series <- Vectorize(function(t,theta,nparams,R)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    series_R <- survival_gen_series_helper(R,nparams)
    1-series_R(t,theta)
}, vectorize.args="t")

#' Survival function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparams[j]} parameters.
#' @param R list of reliability functions
#' @note convert this into a generator over (t,theta)
#' @export
survival_gen_series <- Vectorize(function(t,theta,nparams,R)
{
    stopifnot(length(theta)==sum(nparams))
    series_R <- survival_gen_series_helper(R,nparams)
    function(t,theta) series_R(t,theta)
}, vectorize.args="t")

#' Hazard function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparams[j]} parameters.
#' @param h list of hazard functions
#' @export
hazard_gen_series <- Vectorize(function(t,theta,nparams,h)
{
    stopifnot(length(theta)==sum(nparams))
    series_h <- hazard_gen_series_helper(h,nparams)
    series_h(t,theta)
}, vectorize.args="t")

