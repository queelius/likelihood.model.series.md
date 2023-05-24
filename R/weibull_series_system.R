#' Weibull series
#'
#' This file contains functions related to the Weibull series distribution.
#' Functions include simulation, pdf, cdf, quantile, and other related
#' functions for the Weibull series distribution.
#' 
#' For parameter estimation from masked data, see the file \code{md_weibull_series_mle.R}.
#'
#' @author Alex Towell
#' @name Weibull series
#' @keywords weibull, distribution, series, statistics
NULL

#' Quantile function (inverse of the cdf).
#' By definition, the quantile p * 100% for a strictly monotonically increasing
#' cdf F is the value t that satisfies \code{F(t) - p = 0}.
#' We solve for t using newton's method.
#'
#' @param p vector of probabilities.
#' @param scales vector of weibull scale parameters for weibull lifetime
#'               components
#' @param shapes vector of weibull shape parameters for weibull lifetime
#'               components
#' @param eps stopping condition, default is 1e-3
#' @param t0 initial guess, default is 1
#' @export
qweibull_series <- Vectorize(function(p,scales,shapes,eps=1e-3,t0=1)
{
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    t1 <- NULL
    repeat
    {
        alpha <- 1
        repeat
        {
            t1 <- t0 - alpha * (sum((t0/scales)^shapes) + log(1-p)) /
                sum(shapes*t0^(shapes-1)/scales^shapes)
            if (t1 > 0)
                break
            alpha <- alpha / 2
        }
        if (abs(t1-t0) < eps)
            break
        t0 <- t1
    }
    t1
}, vectorize.args = "p")

#' Sampler for weibull series.
#'
#' NOTE: \code{qweibull_series(p=runif(n),scales,shapes)} is around 6 times slow
#' due to using newton's method to solve \code{F(t) - p = 0}.
#'
#' @param n sample size
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @importFrom stats rweibull
#' @export
rweibull_series <- function(n,scales,shapes)
{
    stopifnot(n > 0)
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    t <- matrix(nrow=n,ncol=m)
    for (j in 1:m)
        t[,j] <- rweibull(n,scale=scales[j],shape=shapes[j])
    apply(t,1,min)
}

#' pdf for weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
dweibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    d <- ifelse(t < 0,
                0,
                sum(shapes/scales*(t/scales)^(shapes-1))*exp(-sum((t/scales)^shapes)))
    ifelse(is.nan(d), 0, d)
}, vectorize.args="t")

#' Survival function for weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
survival_weibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    ifelse(t < 0,
           1,
           exp(-sum((t/scales)^shapes)))
}, vectorize.args="t")

#' Hazard function for weibull series.
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
hazard_weibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    ifelse(t < 0,
           0,
           sum(shapes/scales*(t/scales)^(shapes-1)))
}, vectorize.args="t")


#' The cumulative distribution function for Weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for Weibull component lifetimes
#' @param shapes shape parameters for Weibull component lifetimes
#' @export
pweibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))
    ifelse(t < 0, 0, 1-exp(-sum((t/scales)^shapes)))
}, vectorize.args="t")


#' The cumulative distribution function for Weibull series
#'
#' @param scales scale parameters for Weibull component lifetimes
#' @param shapes shape parameters for Weibull component lifetimes
#' @param g function of system failure time \code{t}, some population
#'          parameter of interest, such as the expected value or the variance.
#'          defaults to expected value, \code{function(t) t}.
#' @importFrom stats integrate
#' @export
param_weibull_series <- function(scales,shapes,g=function(t) t)
{
    integrate(f=function(t) g(t)*dweibull_series(t,scales,shapes),
              lower=0,
              upper=Inf)$value
}
