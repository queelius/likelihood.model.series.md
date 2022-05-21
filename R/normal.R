#' Construct (multivariate or univariate) normal distribution object.
#'
#' @param mu mean
#' @param sigma variance-covariance matrix
#'
#' @export
make_normal <- function(
    mu,
    sigma=diag(length(mu)))
{
    structure(list(
        mu=mu,
        sigma=sigma),
        class=c("normal","dist"))
}

#' Retrieve the variance-covariance matrix (or scalar)
#' of a \code{normal} object.
#'
#' @param object The \code{normal} object to retrieve the variance-covariance matrix from
#' @importFrom stats vcov
#' @export
vcov.normal <- function(object, ...)
{
    object$sigma
}

#' Method for obtaining the parameters of a \code{normal} object.
#'
#' @param x The object to obtain the parameters of
#' @importFrom algebraic.mle params
#' @export
params.normal <- function(x, ...)
{
    c(x$mu,x$sigma)
}

#' Method for sampling from a \code{normal} object.
#'
#' @param x The object to sample from
#' @importFrom mvtnorm rmvnorm
#' @importFrom algebraic.mle sampler
#'
#' @export
sampler.normal <- function(x, ...)
{
    function(n=1)
    {
        mvtnorm::rmvnorm(n,x$mu,x$sigma)
    }
}
