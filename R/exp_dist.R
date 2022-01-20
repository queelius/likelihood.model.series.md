#' Construct exponential series object.
#'
#' @param rate failure rates
#'
#' @export
make_exp_series <- function(rate)
{
    structure(list(
        theta=unlist(rate),
        num_nodes=length(rate)),
        class=c("exp_series","series","exp_dist","distribution"))
}

#' Construct exponential distribution object.
#'
#' @param rate failure rate
#'
#' @export
make_exp_dist <- function(rate)
{
    structure(list(
        theta=unlist(rate),
        num_nodes=length(rate)),
        class=c("exp_dist","random_variable"))
}

#' Method for obtaining the variance-covariance of a \code{exp_series} object.
#'
#' @param exp_dist The object to obtain the variance of
#'
#' @export
vcov.exp_series <- function(exp_series)
{
    diag(1/exp_series$theta)^2
}

#' Method for obtaining the variance of a \code{exp_dist} object.
#'
#' @param exp_dist The object to obtain the variance of
#'
#' @export
vcov.exp_dist <- function(exp_dist)
{
    1/exp_dist$theta^2
}

#' Method for obtaining the parameters of
#' a \code{series} distribution object.
#'
#' @param x The object to obtain the parameters of
#'
#' @export
params.series <- function(x, ...)
{
    x$theta
}

#' Method to obtain the hazard function of
#' an \code{exp_series} object.
#'
#' @param exp_series The object to obtain the hazard function of
#'
#' @export
hazard.exp_dist <- function(exp_dist, ...)
{
    theta <- params(exp_dist)
    function(t,...)
    {
        ifelse(t <= 0,0,sum(theta))
    }
}

#' Method to obtain the pdf of an \code{exp_dist} object.
#'
#' Note that since \code{exp_series} is also exponentially distributed,
#' this works for that too.
#'
#' @param exp_dist The object to obtain the pdf of
#'
#' @export
pdf.exp_dist <- function(exp_dist, ...)
{
    theta <- params(exp_dist)
    function(t,...)
    {
        ifelse(t <= 0,0,sum(theta))
    }
}

#' Method to sample from an \code{exp_dist} object.
#'
#' @param x The object to sample from.
#'
#' @export
sampler.exp_dist <- function(x,...)
{
    function(n=1)
    {
        stats::rexp(params(x))
    }
}
