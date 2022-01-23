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
#' @export
md_fisher_scoring <- function(theta0,info,score,eps=1e-5,max_iterations=250L)
{
    n <- 1L
    repeat
    {
        I <- info(theta0)
        sigma <- ginv(I)
        s <- score(theta0)
        theta1 <- theta0 + sigma %*% s
        if (n == max_iterations || max(abs(theta1-theta0)) < eps)
            return(list(theta.hat=theta1,sigma=sigma,score=s,iterations=n))
        theta0 <- theta1
        n <- n + 1L
    }
}
