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
rexp_system_data = function(n,theta,phi=which.min)
{
    stopifnot(is.numeric(theta))
    data <- system_data(matrix(
        stats::rexp(n*length(theta),rate=theta),ncol=m,byrow=T),
        phi)
    attr(data,"node_lifetime") <- c("exponential")
    data
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
#' @export
dexp_series_failed_node_pdf = function(theta,t)
{
    stopifnot(is.numeric(theta))
    new_parametric_pdf(function(k)
    {
        ifelse(is.wholenumber(k) & k > 0 & k <= length(theta),
               theta[k] / sum(theta), 0)
    },length(theta),1L)
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
loglike_exp_series_masked_data_m0 <- function(masked_data)
{
    stopifnot(!is.null(masked_data$c))
    stopifnot(!is.null(masked_data$s))
    stopifnot("num_nodes" %in% names(attributes(masked_data)))

    m <- masked_data$num_nodes
    new_parametric_func(function(theta)
    {
        theta <- as.matrix(theta,nrow=m,ncol=1)
        sum(log(apply(
            masked_data$c,1,function(r) { sum(theta[r]) }))) -
            sum(masked_data$s) * sum(theta)
    },m,1L)
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
like_exp_series_masked_data_m0 <- function(masked_data)
{
    m <- masked_data$num_nodes
    ll <- loglike_exp_series_masked_data_m0(masked_data)
    new_parametric_func(function(theta) exp(ll(theta)),m,1)
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






mse = function(theta,n) { length(theta)*sum(theta)/n }


info <- function(l,n)
{
    l1 = l[1]
    l2 = l[2]
    l3 = l[3]

    A = matrix(c(1/(l1+l2) + 1/(l1+l3), 1/(l1+l2), 1/(l1+l3),
                 1/(l1+l2), 1/(l1+l2) + 1/(l2+l3), 1/(l2+l3),
                 1/(l1+l3), 1/(l2+l3), 1/(l1+l3) + 1/(l2+l3)),byrow=T,nrow=3)
    n/(2.0*(l1+l2+l3))*A
}


loglike <- function(masked_data)
{
    m <- masked_data$num_nodes
    function(theta)
    {
        sum(log(apply(
            masked_data$c,1,function(r) { sum(theta[r]) }))) -
            sum(masked_data$s) * sum(theta)
    }
}

rmasked_data <- function(n,theta,w)
{
    m <- length(theta)
    t <- matrix(stats::rexp(n*m,rate=theta),ncol=m,byrow=T)
    k <- integer(length=n)
    s <- numeric(length=n)
    for (i in 1:n)
    {
        k[i] <- which.min(t[i,])
        s[i] <- t[i,k[i]]
    }

    c <- matrix(nrow = n, ncol = m, F)

    if (is.null(w)) {
        w <- rep(m, n)
    }
    w <- c(w, rep(w[length(w)], n - length(w)))

    for (i in 1:n)
    {
        c[i, k[i]] <- T
        c[i, sample((1:m)[-k[i]], size = w[i] - 1, replace = F)] <- T
    }
    data <- data.frame(s=s,k=k,t=t,w=w)
    data$c <- c
    data
}





n=1000
theta.star=c(1,3,5)



data=rmasked_data(n,theta.star,2)
ll=loglike(data)
theta.mle=optim(c(1,1,1),function(theta) { -ll(theta) })$par
theta.mle

H.mle=hessian(func=function(theta) { -ll(theta) },x=theta.mle)
#I.mle=info(theta.mle,n)

cov.mle=inv(H.mle)
cov.mle

cov.star=inv(info(theta.star,n))
cov.star

#delta=theta.mle - theta.star

mse.mle = sum(diag(cov.mle))
mse.a = mse(theta.mle,n)
mse.star = mse(theta.star,n)

mse.mle
mse.star
mse.a



#### sufficient statistics

d = rmasked_data(10,theta.star,2)
loglike_ss = function(masked_data)
{
    n = nrow(masked_data)
    m = ncol(masked_data$c)
    freq = d %>% group_by(c) %>% count

    #t(apply(freq$c,1,function(r) { (1:m)[r]}))

}

loglike_ss(data)
