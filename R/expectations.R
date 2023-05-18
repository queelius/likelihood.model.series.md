#' Expectation operator with respect to pdf `f` of function `g`,
#' i.e., `E(g(T))`, assuming `T ~ f(t)`.
#'
#' @param f pdf
#' @param g characteristic function of interest, defaults to T
#' @param lower lower limit of support of `f`, defaults to 0
#' @param upper upper limit of support of `f`, defaults to infinity
#' @importFrom cubature adaptIntegrate
#' @export
expectation <- function(f,g=function(t) t, lower=0,upper=Inf)
{
    adaptIntegrate(function(t) f(t)*g(t),lowerLimit=lower,upperLimit=upper)$integral
}

#' Monte carlo expectation of a random vector `T ~ sampler(1)` applied to a
#' function `g` where `g` should have a domain that is a superset of
#' the random vector's support. That is, we estimate `E(g(T))` with the
#' mean of `g(T[1]),...,g(T[n])`. We see that variance is
#' `1/n var(g(T[1]))`, which is inversely proportional to the sample size
#' `n`.
#'
#' @param sampler sampler, generates random vectors for a desired distribution
#' @param g characteristic function of interest, defaults to identity
#' @param n number of sample points, larger samples reduce variance
#' @export
mc_expectation <- function(sampler,g=function(t) t,n=10000)
{
    y <- g(sampler(n))
    ifelse(is.numeric(y), mean(y), colMeans(y))
}

