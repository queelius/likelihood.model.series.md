#' Compute the covariance matrix from the given masked data estimate.
#'
#' Sampling distribution of the MLE is a multivariate normal with mean
#' given by the true parameter value and, asymptotically, a covariance
#' given by the inverse of the Fisher information matrix.
#' @export
vcov.md_estimate <- function(estimate)
{
    matlib::inv(estimate$info)
}

#' Method to obtain the confidence intervals of
#' the parameter values of a masked data
#' estimator, \code{md_estimate}.
#'
#' @param level Confidence level, defaults to 0.95 (alpha=.05)
#'
#' @export
confint.md_estimate <- function(estimate, level=0.95)
{
    V <- vcov(estimate)
    theta.hat <- point(estimate)
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
#' @export
point <- function(x, ...)
{
    UseMethod("point",x)
}

#' Method to obtain the point estimate of
#' a masked data estimator, \code{md_estimate}.
#'
#' @export
point.md_estimate <- function(estimate)
{
    estimate$theta.hat
}
