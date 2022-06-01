#' Maximum likelihood estimator of the parameters of a series
#' system with nodes that have exponentially distributed
#' lifetimes given a sample of masked data according to
#' candidate model m0.
#'
#' @param md masked data
#' @param loglike log-likelihood function
#' @param theta0 initial guess for MLE
#' @param eps stopping condition, defaults to \code{1e-05}
#' @param max_iter stop if iterations reaches max_iter, default is 0 (infinity)
#' @return numeric MLE estimate
#' @importFrom algebraic.mle mle_gradient_ascent
#' @export
md_grad_mle_series_C1_C2_C3 <- function(
        md,
        loglike,
        theta0,
        eps=1e-5,
        max_iter=0L)
{
    mle_gradient_ascent(
        l=loglike(md),
        theta0=theta0,
        sup=function(theta) all(theta > 0),
        stop_cond=function(theta1,theta0) abs(max(theta1-theta0)) < eps,
        max_iter=max_iter)
}
