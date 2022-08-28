#' Expectation operator with respect to pdf \code{f} of function \code{g},
#' i.e., \code{E(g(T))}, assuming \code{T ~ f(t)}.
#'
#' @param f pdf
#' @param g characteristic function of interest, defaults to T
#' @param lower lower limit of support of \code{f}, defaults to 0
#' @param upper upper limit of support of \code{f}, defaults to infinity
#' @importFrom cubature adaptIntegrate
#' @export
expectation <- function(f,g=function(t) t, lower=0,upper=Inf)
{
    adaptIntegrate(function(t) f(t)*g(t),lowerLimit=lower,upperLimit=upper)$integral
}

#' Monte carlo expectation of a random vector \code{T ~ sampler(1)} applied to a
#' function \code{g} where \code{g} should have a domain that is a superset of
#' the random vector's support. That is, we estimate \code{E(g(T))} with the
#' mean of \code{g(T[1]),...,g(T[n])}. We see that variance is
#' \code{1/n var(g(T[1]))}, which is inversely proportional to the sample size
#' \code{n}.
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

# Retrieve the function arguments.
md_func_args <- function(...)
{
    call <- evalq(match.call(expand.dots = F), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))

    for(i in setdiff(names(formals), names(call)))
        call[i] <- list( formals[[i]] )

    match.call(sys.function(sys.parent()), call)
}

# mean absolute deviation
mad <- function(x)
{
    u <- mean(x)
    n <- length(x)
    sum(abs(x-u)) / n
}

# linearly rescale numbers in x each in the range [x_low,x_high] to an output
# in the range [out_low, out_high].
rescale <- function(x,x_low=0,x_high=1,y_low=0,y_high=1)
{
    y_low + ((y_high - y_low) / (x_high - x_low)) * (x - x_low)
}

survival_general_series_helper <- function(R,nparams)
{
    m <- length(R)
    function(t,theta)
    {
        index.0 <- 1
        index.1 <- nparams[1]
        p <- 0
        for (j in 1:m)
        {
            p <- p + R[[j]](t,theta[index.0:index.1])
            index.0 <- index.1 + 1
            index.1 <- index.1 + nparams[j]
        }
        p
    }
}

hazard_general_series_helper <- function(h,nparams)
{
    m <- length(h)
    function(t,theta)
    {
        index.0 <- 1
        index.1 <- nparams[1]
        v <- 0
        for (j in 1:m)
        {
            v <- v + h[[j]](t,theta[index.0:index.1])
            index.0 <- index.1 + 1
            index.1 <- index.1 + nparams[j]
        }
        v
    }
}

cumulative_hazard <- function(h,t)
{
    adaptIntegrate(h, lowerLimit=0, upperLimit=t)$integral
}

generate_survival_from_hazard <- function(h)
{
    function(t) exp(-cumulative_hazard(h,t))
}

