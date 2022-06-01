#' Masked data approximately satisfies the following set of conditions:
#' C1: \code{Pr(K[i] in C[i]) = 1}
#' C2: \code{Pr(C[i]=c[i] | K[i]=j, T[i]=t[i]) = Pr(C[i]=c[i] | K[i]=j', T[i]=t[i])}
#'     for any \code{j,j' in c[i]}.
#' C3: masking probabilities are independent of theta
#'
#' @param md masked data
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of survival functions
#' @returns a log-likelihood function with respect to theta given \code{md}
#' @importFrom md.tools md_decode_matrix
#' @export
md_loglike_general_series_C1_C2_C3 <- function(md,nparams,h,R)
{
    m <- length(h)
    stopifnot(m > 0)
    stopifnot(m==length(R))
    series_h <- hazard_general_series_helper(h,nparams)
    series_R <- survival_general_series_helper(R,nparams)

    C <- md_decode_matrix(md,"x")
    stopifnot(ncol(C)==m)
    n <- nrow(md)
    stopifnot(n > 0)

    right_censoring <- "delta" %in% colnames(md)
    stopifnot(!right_censoring || "s" %in% colnames(md))
    t <- ifelse(right_censoring, md$s, md$t)

    log.R <- function(t,theta)
    {
        sum <- 0
        for (Rj in R)
            sum <- sum + log(Rj(t,theta))
        sum
    }

    function(theta)
    {
        stopifnot(length(theta)==sum(nparams))
        res <- 0
        for (i in 1:n)
        {
            if (right_censoring && md$delta[i])
            {
                res <- res + log.R(t[i],theta)
            }
            else #if (!right_censoring || (right_censoring && !md$delta[i]))
            {
                haz <- 0
                for (j in (1:m)[C[i,]])
                    haz <- haz + h[[j]](t[i],theta)
                res <- res + log(haz) * log.R(t[i],theta)
            }
        }
        res
    }
}


#' Quantile function (inverse of the cdf) for a general series system.
#' By definition, the quantile \code{p * 100%} is the value \code{t} that
#' satisfies \code{F(t) - p = 0}. We solve for \code{t} using Newton's method.
#'
#' @param p vector of probabilities.
#' @param theta parameter vector to evaluate h's and R's at
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of reliability functions
#' @param eps stopping condition, default is 1e-3
#' @param t0 initial guess, default is 1
#' @export
qgeneral_series <- Vectorize(function(p,theta,nparams,h,R,eps=1e-3,t0=1)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    stopifnot(eps > 0)
    stopifnot(all(p > 0))

    series_h <- hazard_general_series_helper(h,nparams)
    series_R <- survival_general_series_helper(R,nparams)
    t1 <- NULL
    repeat
    {
        alpha <- 1
        repeat
        {
            t1 <- t0 - alpha * (1-series_R(t,theta))/(series_h(t,theta)*series_R(t,theta))
            if (t1 > 0)
                break
            alpha <- alpha / 2
        }
        if (abs(t2-t1) < eps)
            break
        t1 <- t2
    }
    t2
}, vectorize.args = "p")

#' Sample from a general series system of \code{m} components whose hazard
#' and reliability functions are respectively given by \code{h} and \code{R}.
#'
#' @param n sample size
#' @param theta parameter vector to evaluate h's and R's at
#' @param h list of hazard functions
#' @param R list of reliability functions
#' @export
rgeneral_series <- function(n,theta,h,R)
{
    qgeneral_series(runif(n),theta,h,R)
}

#' pdf for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector to evaluate h's and R's at
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of reliability functions
#' @export
dgeneral_series <- Vectorize(function(t,theta,nparams,h,R)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    series_h <- hazard_general_series_helper(h,nparams)
    series_R <- survival_general_series_helper(R,nparams)
    series_h(t,theta) / series_R(t,theta)
}, vectorize.args="t")

#' cdf for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector to evaluate R's at
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param R list of reliability functions
#' @export
pgeneral_series <- Vectorize(function(t,theta,nparams,R)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    series_R <- survival_general_series_helper(R,nparams)
    1-series_R(t,theta)
}, vectorize.args="t")

#' Survival function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector to evaluate R's
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparams[j]} parameters.
#' @param R list of reliability functions
#' @export
survival_general_series <- Vectorize(function(t,theta,nparams,R)
{
    stopifnot(length(theta)==sum(nparams))
    series_R <- survival_general_series_helper(R,nparams)
    series_R(t,theta)
}, vectorize.args="t")

#' Hazard function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector to evaluate h's at
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparams[j]} parameters.
#' @param h list of hazard functions
#' @export
hazard_general_series <- Vectorize(function(t,theta,nparams,h)
{
    stopifnot(length(theta)==sum(nparams))
    series_h <- hazard_general_series_helper(h,nparams)
    series_h(t,theta)
}, vectorize.args="t")

