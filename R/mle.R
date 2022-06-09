#' Maximum likelihood estimator of the parameters of a series
#' system given a sample of masked data according to candidate sets
#' that are generated in a way that is compatible with conditions C1, C2, and
#' C3 so that we do not need to know or estimate the distribution of candidate
#' sets.
#'
#' @param loglike log-likelihood function
#' @param theta0 initial guess for MLE
#' @param theta.restart a function \code{() -> R^m} for generating new theta0 for restarts, defaults to randomly generating values
#' @param sup support function, defaults to all components of \code{theta} being positive non-zero
#' @param eps stopping condition, defaults to \code{1e-05}
#' @param restarts number of restarts to do, defaults to 100
#' @param max_iter stop if iterations reaches max_iter, default is 0 (infinity)
#' @param eta learning rate, defaults to 1
#' @param r backtracking line search parameter
#' @return numeric MLE estimate
#' @importFrom algebraic.mle mle_gradient_ascent
#' @importFrom dplyr near
#' @export
md_mle_series_C1_C2_C3 <- function(
        loglike,
        theta0,
        theta.restart=NULL,
        sup=function(theta) all(theta > 0),
        eps=1e-5,
        restarts=100L,
        max_iter=250L,
        eta=1,
        r=0.5)
{
    m <- length(theta0)
    if (is.null(theta.restart))
        theta.restart <- function() rexp(n=m,rate=.2)
    max_loglik <- -Inf
    mles <- list()
    for (i in 1:restarts)
    {
        mle <- mle_gradient_ascent(
            l=loglike,
            theta0=theta0,
            sup=sup,
            stop_cond=function(theta1,theta0) abs(max(theta1-theta0)) < eps,
            max_iter=max_iter)

        loglik <- loglike(mle)
        if (near(loglik,max_loglik))
        {
            is_near <- T
            for (x in mles)
            {
                if (!all(near(point(x),point(mle))))
                {
                    is_near <- F
                    break
                }
            }
            if (is_near)
            {
                n <- length(mles)
                max_loglik <- ifelse(n==0,loglik,max_loglik+(loglik-max_loglik)/n)
                mles <- append(mles,mle)
            }
        }
        else if (loglik > max_loglik)
        {
            mles <- list(mle)
            max_loglik <- loglik
        }

        if (i != restarts)
        {
            repeat
            {
                theta0 <- theta.restart()
                if (sup(theta0))
                    break
            }
        }
    }
}



# beta: shape factor, controls the type of failure of the element (infant mortality, wear-out, or random).
# beta = k
# eta: scale factor, representing the time when 63.2 % of the total population is failed.
# eta = lambda
