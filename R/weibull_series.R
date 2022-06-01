
#' Masked data approximately satisfies the following set of conditions:
#' C1: \code{Pr(K[i] in C[i]) = 1}
#' C2: \code{Pr(C[i]=c[i] | K[i]=j, T[i]=t[i]) = Pr(C[i]=c[i] | K[i]=j', T[i]=t[i])
#'     for any j, j' in c[i]}.
#' C3: masking probabilities are independent of theta
#'
#' @param md masked data
#' @returns a log-likelihood function with respect to theta given md
#' @importFrom md.tools md_decode_matrix
#' @export
md_loglike_weibull_series_C1_C2_C3 <- function(md)
{
    C <- md_decode_matrix(md,"x")
    m <- ncol(C)
    n <- nrow(md)
    stopifnot(m > 0)
    stopifnot(n > 0)

    function(theta)
    {
        # theta should be a parameter vector of length 2*m
        # todo: make first m values of `theta` be for scale parameters
        #            second m values (m+1,...,2m) be for shape parameters
        scales <- theta[(0:(m-1)*2)+1]
        shapes <- theta[(1:m)*2]
        s <- 0

        for (i in 1:n)
        {
            acc <- 0
            for (j in 1:m)
                acc <- acc + (md$s[i]/scales[j])^shapes[j]
            s <- s + log(acc)

            if (!md$delta[i])
            {
                acc <- 0
                c <- (1:m)[C[i,]]
                for (j in c)
                    acc <- acc + shapes[j]/scales[j]*(md$s[i]/scales[j])^(shapes[j]-1)
                s <- s + log(acc)
            }
        }
        s
    }
}


#' Quantile function (inverse of the cdf).
#' By definition, the quantile p * 100% for a strictly monotonically increasing
#' cdf F is the value t that satisfies \code{F(t) - p = 0}.
#' We solve for t using newton's method.
#'
#' @param p vector of probabilities.
#' @param scales vector of weibull scale parameters for weibull lifetime
#'               components
#' @param shapes vector of weibull shape parameters for weibull lifetime
#'               components
#' @param eps stopping condition, default is 1e-3
#' @param t0 initial guess, default is 1
#' @export
qweibull_series <- Vectorize(function(p,scales,shapes,eps=1e-3,t0=1)
{
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    t1 <- NULL
    repeat
    {
        alpha <- 1
        repeat
        {
            t1 <- t0 - alpha * (sum((t0/scales)^shapes) + log(1-p)) /
                sum(shapes*t0^(shapes-1)/scales^shapes)
            if (t1 > 0)
                break
            alpha <- alpha / 2
        }
        if (abs(t1-t0) < eps)
            break
        t0 <- t1
    }
    t1
}, vectorize.args = "p")

#' Sampler for weibull series.
#'
#' NOTE: \code{qweibull_series(p=runif(n),scales,shapes)} is around 6 times slow
#' due to using newton's method to solve \code{F(t) - p = 0}.
#'
#' @param n sample size
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @importFrom stats rweibull
#' @export
rweibull_series <- function(n,scales,shapes)
{
    stopifnot(n > 0)
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    t <- matrix(nrow=n,ncol=m)
    for (j in 1:m)
        t[,j] <- rweibull(n,scale=scales[j],shape=shapes[j])
    apply(t,1,min)
}

#' pdf for weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
dweibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    ifelse(t < 0,
           0,
           sum(shapes/scales*(t/scales)^(shapes-1))*exp(-sum((t/scales)^shapes)))
}, vectorize.args="t")

#' Survival function for weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
survival_weibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    ifelse(t < 0,
           1,
           exp(-sum((t/scales)^shapes)))
}, vectorize.args="t")

#' Hazard function for weibull series.
#'
#' @param t series system lifetime
#' @param scales scale parameters for weibull component lifetimes
#' @param shapes shape parameters for weibull component lifetimes
#' @export
hazard_weibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))

    ifelse(t < 0,
           0,
           sum(shapes/scales*(t/scales)^(shapes-1)))
}, vectorize.args="t")


#' The cumulative distribution function for Weibull series
#'
#' @param t series system lifetime
#' @param scales scale parameters for Weibull component lifetimes
#' @param shapes shape parameters for Weibull component lifetimes
#' @export
pweibull_series <- Vectorize(function(t,scales,shapes)
{
    m <- length(scales)
    stopifnot(m==length(shapes))
    stopifnot(all(shapes > 0))
    stopifnot(all(scales > 0))
    ifelse(t < 0, 0, 1-exp(-sum((t/scales)^shapes)))
}, vectorize.args="t")




