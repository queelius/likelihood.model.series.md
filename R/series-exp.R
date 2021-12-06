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

series.exp.mle = function(masked.data)
{
    # closed-form solution exists, see paper
    1
}

series.exp.cov = function(masked.data)
{
    # first, estimate mle using series.exp.mle
    # then, either use the expected info or the
    # observed info. observed info is just
    #   -hessian(loglike.series.exp(rate.mle,masked.data))
    # a nice feature of the observed info is we
    # allow the joint distribution of C and alpha
    # to be a function of the masked.data.
    1
}
