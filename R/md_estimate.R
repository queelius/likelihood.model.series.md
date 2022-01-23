#' Compute the covariance matrix from the given masked data estimate.
#'
#' Sampling distribution of the MLE is a multivariate normal with mean
#' given by the true parameter value and, asymptotically, a covariance
#' given by the inverse of the Fisher information matrix.
#'
#' @param object The variance-covariance matrix of the estimator to obtain
#'
#' @export
vcov.md_estimate <- function(object, ...)
{
    object$sigma
}

#' Method to obtain the confidence intervals of the parameter values of a
#' masked data estimator, \code{md_estimate}.
#'
#' @param object The \code{md_estimate} object to compute the confidence intervals for
#' @param parm Unused
#' @param level Confidence level, defaults to 0.95 (alpha=.05)
#'
#' @export
confint.md_estimate <- function(object, parm=NULL, level=0.95, ...)
{
    V <- vcov(object)
    theta.hat <- point(object)
    p <- length(theta.hat)
    q <- stats::qnorm(level)

    ci <- matrix(rep(NA,p*2),p,2)
    colnames(ci) <- c(paste((1-level)/2*100,"%",sep=""),
                      paste((1-(1-level)/2)*100,"%",sep=""))
    for (j in 1:p)
    {
        ci[j,1] <- theta.hat[j] - q * sqrt(V[1,1])
        ci[j,2] <- theta.hat[j] + q * sqrt(V[1,1])
    }
    ci
}

#' Generic method for obtaining the point estimate of
#' an estimator.
#'
#' @param x The object to obtain the point estimate of
#'
#' @export
point <- function(x, ...)
{
    UseMethod("point",x)
}

#' Method to obtain the point estimate of
#' a masked data estimator, \code{md_estimate}.
#'
#' @param x The \code{md_estimate} object to obtain the point estimate of
#'
#' @export
point.md_estimate <- function(x, ...)
{
    x$theta.hat
}

#' Generic method for obtaining the fisher information
#' matrix of an estimator, if supported.
#'
#' @param x The object to obtain the fisher information of
#'
#' @export
fisher_info <- function(x, ...)
{
    UseMethod("fisher_info",x)
}

#' Method to obtain the fisher information matrix of an \code{md_estimate}.
#'
#' @param x The \code{md_estimate} object to obtain the fisher information of
#'
#' @export
fisher_info.md_estimate <- function(x, ...)
{
    x$info
}

#' Method to obtain the sampler for an \code{md_estimate} object.
#'
#' @param x The \code{md_estimate} object to create a sampling procedure from
#' @export
sampler.md_estimate <- function(x, ...)
{
    sampler(
        make_normal(
            point(x),
            vcov(x)))
}
