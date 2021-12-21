#' @brief rseries.masked.exp generates the joint distribution
#'        of a series system with exponentially distributed
#'        components.
#'
#' @param n sample size (each row is an observation)
#' @param rate the j-th component has a failure rate rate_j
#'
#' @return dataframe of n observations, S x K x T. K and T
#'         are latent, and some other covariate or predictor
#'         will be necessary to estimate theta.
#' @export
#'
#' @examples
#' # generate 10 samples
#' masked.data = rseries.exp(n=10,theta=c(1,2,3))
#' # get system failure times only
#' lifetimes = masked.data$s
rseries.masked.exp = function(n,theta)
{
    m = length(theta)
    t = matrix(stats::rexp(n*m,rate=theta),ncol=m,byrow=T)
    df = data.frame(s=apply(t,1,min),
                    k=apply(t,1,which.min))
    df$t = t

    class(df) = c("masked.data",class(df))
    attr(df,nodes) = length(theta)
    attr(df,latent) = c("k","t")
    df
}

#' @brief deseries.node.exp probability mass
#'        function for f(k|t,theta).
#'
#' note that f(k) = integrate f(k,t) dt
#' is the same as f(k,t) for the
#' exponential series case, due to
#' the memoryless property of
#' the exponential distribution.
#'
#' @param k
#' @param t
#' @return probability that the failed
#'         component is k at time t.
#'
#' @export
#'
dseries.node.exp = function(theta,t)
{
    stopifnot(is.numeric(theta))
    function(k)
    {
        ifelse(is.wholenumber(k) & k > 0 & k <= length(theta),
               theta[k] / sum(theta), 0)
    }
}

#' @brief log-likelihood function of masked data for a series system
#'        with exponentially distributed lifetimes.
#'
#' @param masked.data alpha-masked data
#'
#' @return log-likelihood function
#'             l : R^m -> R
#'                 theta |-> log-likelihood(theta | masked.data)
#' @export
llseries.masked.exp <- function(masked.data)
{
    function(theta)
    {
        sum(log(apply(
            masked.data$c,1,function(r) { sum(theta[r]) }))) -
            sum(masked.data$s) * sum(theta)
    }
}

#' likelihood function of masked data
#' for a series system with exponentially distributed lifetimes.
#'
#' @param masked.data masked data
#'
#' @return likelihood function
#'             L : R^m -> R
#'                 theta |-> likelihood(theta | masked.data)

#' @export
lseries.exp <- function(masked.data)
{
    exp(llseries.exp(masked.data))
}

#' maximum likelihood estimator of the
#' parameters of a series system with
#' components that have exponentially
#' distributed lifetimes given a
#' sample of masked data.
#'
#' @param masked.data masked data
#'
#' @return mle of theta given masked data
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
#' @param theta rate parameter
#'
#' @return expected info
#' @export
series.exp.mle.observed.info = function(masked.data)
{
    ll = llseries.exp(masked.data)
    function(theta)
    {
        hessian(func=ll,x=theta)
    }
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
#' @param theta rate parameter
#'
#' @return expected info
#' @export
series.exp.mle.info = function(n,theta)
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
series.exp.mle.cov = function(n,theta)
{
    series.exp.mle.info(n,theta)^(-1)
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
