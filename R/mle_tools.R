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
#' @param gamma step size
#'
#' @return MLE estimate of theta
#' @author Alex Towell
#' @export
md_fisher_scoring <- function(theta0,info,score,eps=1e-5,max_iterations=10000L,gamma=.1)
{
    n <- 1L
    repeat
    {
        theta1 <- theta0 + gamma * matlib::inv(info(theta0)) %*% score(theta0)
        if (n == max_iterations || max(abs(theta1-theta0)) < eps)
            return(list(theta.hat=theta1,iterations=n))
        theta0 <- theta1
        n <- n + 1L
    }
}

#' Fisher scoring algorithm using the log-likelihood (or kernel log-like)
#' function.
#'
#' @param theta0 initial guess of theta with p components
#' @param loglike information matrix function of type R^p -> R^(p x p)
#' @param eps stopping condition
#'
#' @return MLE estimate of theta
#' @export
md_fisher_scoring_loglike <- function(theta0,loglike,eps=1e-5,max_iterations=10000L,gamma=.1)
{
    md_fisher_scoring(
        theta0,
        function(theta) { numDeriv::hessian(function(x) { -loglike(x) },theta) },
        function(theta) { numDeriv::grad(loglike,theta) },
        eps,
        max_iterations,
        gamma)
}





#' Gradient ascent algorithm.
#'
#' @section Algorithm:
#' The algorithm is straightforward. Details here.
#'
#' @param theta0 initial guess of theta with \eqn{p} components
#' @param score score function of type \eqn{R^p -> R^p}
#' @param eps stopping condition
#' @param max_iterations maximum number of iterations
#' @param gamma step size
#' @author Alex Towell
#' @export
md_solver_gradient_ascent <- function(theta0,grad,eps=1e-5,max_iterations=10000L,gamma=.1)
{
    n <- 1L
    repeat
    {
        theta1 <- theta0 - grad(theta0)
        if (n == max_iterations || max(abs(theta1-theta0)) < eps)
            return(list(theta.hat=theta1,iterations=n))
        theta0 <- theta1
        n <- n + 1L
    }
}

