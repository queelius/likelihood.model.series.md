#' Fisher scoring algorithm.
#'
#' @section Algorithm:
#' The algorithm is straightforward. Details here.
#'
#' @param theta0 initial guess of theta with \eqn{p} components
#' @param info information matrix function of type \eqn{R^p -> R^{p \times q}}
#' @param score score function of type \eqn{R^p -> R^p}
#' @param eps stopping condition
#' @param max_iterations maximum number of iterations
#'
#' @return MLE estimate of theta
#' @author Alex Towell
#' @export
md_fisher_scoring <- function(theta0,info,score,eps=1e-5,max_iterations=10000L)
{
    n <- 1L
    repeat
    {
        theta1 <- theta0 + matlib::inv(info(theta0)) %*% score(theta0)
        if (n == max_iterations || max(abs(theta1-theta0)) < eps)
            return(list(theta.hat=theta1,iterations=n, max_iterations=max_iterations==n))
        theta0 <- theta1
        n <- n + 1L
    }
}

#' Fisher scoring algorithm using the log-likelihood function (or more ideally, kernel)
#'
#' @param theta0 initial guess of theta with p components
#' @param loglike information matrix function of type R^p -> R^(p x p)
#' @param eps stopping condition
#'
#' @return MLE estimate of theta
#' @export
md_fisher_scoring_loglike <- function(theta0,loglike,eps=1e-5)
{
    md_fisher_scoring(
        theta0,
        function(theta) { numDeriv::hessian(function(x) { -loglike(x) },theta) },
        function(theta) { numDeriv::grad(loglike,theta) })
}

#' Compute the covariance matrix from the information matrix for
#' an MLE parameter.
#'
#' Sampling distribution of the MLE is
#' a multivariate normal with mean
#' given by the true parameter value
#' and, asymptotically, a covariance
#' given by the inverse of the Fisher
#' information matrix.
#'
#' @param info information matrix
#' @export
md_mle_covariance_fisher <- function(info)
{
    function(theta) { numDeriv::inv(info(theta)) }
}


#' Compute the covariance matrix from a sample
#' of MLEs.
#'
#' The asymptotic sampling distribution of the MLE
#' is a multivariate normal with mean
#' given by the true parameter value
#' and a covariance given by the inverse of its
#' Fisher information matrix.
#'
#' However, if we are not at the asymptotic
#' limit, we may estimate the sampling distribution
#' if we have ...
#'
#' @param mles vector of MLEs
#' @export
md_mle_covariance <- function(mles)
{
    NA
}
