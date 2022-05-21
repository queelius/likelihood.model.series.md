#' Construct exponential distribution object.
#'
#' @param rate failure rate
#'
#' @export
make_exp_dist <- function(rate)
{
    structure(list(
        theta=unlist(rate),
        num_comp=length(rate)),
        class=c("exp_dist","dist"))
}


#' Method for obtaining the variance of a \code{exp_dist} object.
#'
#' @param object The \code{exp_dist} object to obtain the variance of
#' @importFrom stats vcov
#' @export
vcov.exp_dist <- function(object, ...)
{
    1/object$theta^2
}

#' Method for obtaining the parameters of
#' a \code{exp_dist} distribution object.
#'
#' @param x The \code{exp_dist} object to obtain the parameters of
#' @importFrom algebraic.mle params
#' @export
params.exp_dist <- function(x, ...)
{
    x$theta
}

#' Method to obtain the hazard function of
#' an \code{exp_dist} object.
#'
#' @param x The \code{exp_dist} object to obtain the hazard function of
#'
#' @export
hazard.exp_dist <- function(x, ...)
{
    theta <- params(x)
    function(t,...)
        ifelse(t <= 0,0,theta)
}



#' Method to obtain the pdf of an \code{exp_dist} object.
#'
#' @param x The object to obtain the pdf of
#'
#' @export
pdf.exp_dist <- function(x, ...)
{
    theta <- params(x)
    function(t,...)
    {
        ifelse(t <= 0,0,theta)
    }
}

#' Method to sample from an \code{exp_dist} object.
#'
#' @param x The \code{exp_dist} object to sample from.
#' @importFrom algebraic.mle sampler
#' @export
sampler.exp_dist <- function(x,...)
{
    function(n=1,...)
    {
        stats::rexp(n,params(x),...)
    }
}


#' Method to obtain the fisher information of the \code{x} object that
#' represents an \code{exp_dist} object.
#'
#' @importFrom algebraic.mle fisher_info
#' @export
fisher_info.exp_dist <- function(x,...)
{
    1/params(x)^2
}

