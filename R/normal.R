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
        class=c("normal","random_variable"))
}

#' Retrieve the variance-covariance matrix (or scalar)
#' of a \code{normal} object.
#'
#' @param object The \code{normal} object to retrieve the variance-covariance matrix from
#'
#' @export
vcov.normal <- function(object, ...)
{
    object$sigma
}

#' Method for obtaining the parameters of a \code{normal} object.
#'
#' @param x The object to obtain the parameters of
#'
#' @export
params.normal <- function(x, ...)
{
    c(x$mu,x$sigma)
}

#' Method for sampling from a \code{normal} object.
#'
#' @param x The object to sample from
#'
#' @export
sampler.normal <- function(x, ...)
{
    function(n=1)
    {
        mvtnorm::rmvnorm(n,x$mu,x$sigma)
    }
}
