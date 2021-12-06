#' rseries.exp
#'
#' @param n sample size (each row is an observation)
#' @param rate the j-th component has a failure rate rate_j
#'
#' @return matrix of n x length(rates) component lifetimes
#' @export
#'
#' @examples
#' # generate 10 samples (10 x 3 matrix of component lifetimes)
#' t = rseries.exp(n=10,rate=c(1,2,3))
rseries.exp = function(n,rate)
{
    m = length(rate)
    t = matrix(stats::rexp(n*m,rate=rate),ncol=m,byrow=T)
    K = apply(t,1,which.min)
    S = apply(t,1,min)
    data.frame(S=S,K=K,t=t)
}

#' log-likelihood function of masked data
#' for a series system with exponentially distributed lifetimes.
#'
#' @param rate rate parameter
#' @param masked.data masked data
#'
#' @return log-likelihood of masked data wrt rate parameter
#' @export
loglike.series.exp = function(rate,masked.data)
{
    m = ncol(masked.data)-2
    Cs = masked.data[,3:m]
    n = nrow(Cs)
    s = 0
    for (i in 1:n)
        s = s + log(sum(rate[Cs[i,]]))
    s - sum(masked.data$S) * sum(rate)
}

#' maximum likelihood estimator of the
#' parameters of a series system with
#' components that have exponentially
#' distributed lifetimes given a
#' sample of masked data.
#'
#' @param masked.data masked data
#'
#' @return mle of rate parameter
#' @export
series.exp.mle = function(masked.data)
{
    # closed-form solution exists, see paper
    1
}

#' information matrix of the rate parameter
#' of the series system with exponentially
#' distributed component lifetimes. this is
#' the expected info. if observed info is
#' desired, then just evaluate the hessian
#' of the log-likelihood function at the
#' mle.
#'
#' @param n masked data sample size
#' @param rate true rate (or mle)
#'
#' @return expected info
#' @export
series.exp.info = function(n,rate)
{
    1
}

#' sampling distribution of the mle is
#' a multivariate normal with mean
#' given by the true rate parameterfor a series system
#' with exponentially distributed component
#' lifetimes, given a sample of masked data.
#'
#' @param n sample size
#' @param rate true rate parameter (or mle)
#'
#' @return multivariate normal of the mle's sampling distribution
#' @export
series.exp.mle.cov = function(n,rate)
{
    series.exp.info(n,rate)^(-1)
}

#' bootstrap of the mle's sampling distribution
#'
#' @param masked.data sample of masked data
#' @param r replicates
#'
#' @return mle and bootstrap estimate of its covariance matrix
#' @export
rseries.exp.mle.cov.bootstrap = function(masked.data,r)
{
    n = nrow(masked.data)
    rate.mle = series.exp.mle(masked.data)
    rate.bs = NULL
    for (i in 1:r)
    {
        d = rmasked.data(n,rate.mle)
        rate.bs = rbind(rate.bs,series.exp.mle(d))
    }
    list(mle=rate.mle,Sigma=cov(rate.bs))
}
