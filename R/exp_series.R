#' @brief rexp.series.m0 generates the joint distribution
#'        of a series system with exponentially distributed
#'        components according to model m0
#'
#' @param n sample size (each row is an observation)
#' @param theta the j-th component has a failure rate theta[j]
#'
#' @return data frame of n observations, (s,k,t1,...,tm,c1,...,cm)
#'         where k's and t's and c's are a covariate or predictor
#'         of theta.
#' @export
#'
#' @examples
#' # generate 10 samples
#' masked.data = rseries.exp(n=10,theta=c(1,2,3))
#' # get system failure times only
#' lifetimes = masked.data$s
# preconditions: n >= 0,
#                theta[j] > 0 for all j
#                w[j] > 0 for all j
rexp.series.m0 <- function(n,theta,w)
{
    m <- length(theta)

    t <- as_tibble(matrix(stats::rexp(n*m,rate=theta), ncol=m, byrow=T))
    names(t) <- paste("t",1:m,sep=".")

    md <- tibble(s = apply(t,1,min),
                 k = apply(t,1,which.min),
                 w = w)
    md <- rcandidates.m0(md)
    md <- bind_cols(md,t)

    # set up attribute properties
    attr.nodes <- list()
    for (j in 1:m)
        attr.nodes[[j]] <- list("family" = "exponential",
                                "index"  = theta[j])

    attr(md,"sim") <- list("theta" = theta,
                           "nodes" = attr.nodes,
                           "data"  = list("model" = "m0"))
    attr(md,"masked") <- c("k",paste("t",1:m,sep="."))
}

#' @brief log-likelihood function of masked data for a series system
#'        with exponentially distributed lifetimes.
#'
#' @param md masked data for candidate model m0
#'
#' @return log-likelihood function
#'             l : R^m -> R
#'                 theta |-> log-likelihood(theta | masked.data)
#' @export
loglike.exp.series.m0 <- function(md)
{
    m <- md.nnodes(md)
    C <- md.candidates.matrix(md)

    function(theta)
    {
        theta <- as.matrix(theta,nrow=m,ncol=1)
        sum(log(apply(
            C,1,function(r) { sum(theta[r]) }))) -
            sum(md$s) * sum(theta)
    }
}

# loglike.exp.series.m0 <- function(md)
# {
#     library(dplyr)
#     library(matlib)
#
#     md$C <- md.candidates.matrix(md)
#     freq = md %>% group_by(C) %>% count
#     sum.t = -sum(md$s)
#     n = nrow(freq)
#     function(theta)
#     {
#         f <- sum.t * sum(theta)
#         for (i in 1:n)
#         {
#             f <- f + freq$n[i] * log(sum(theta[freq$C[i,]]))
#         }
#         f
#     }
# }

score.series.m0 <- function(md,p)
{
    md$C <- md.candidates.matrix(md)
    freq = md %>% group_by(C) %>% count
    sum.t = -sum(md$s)
    p = ncol(md$C)

    function(theta)
    {
        v <- matrix(rep(sum.t,p),nrow=p)
        for (j in 1:p)
        {
            f = freq %>% filter(C[[j]] == T)
            for (i in 1:nrow(f))
            {
                S <- sum(theta[f$C[i,]])
                v[j,1] <- v[j,1] + f$n[i] / S
            }
        }
        v
    }
}

info.series.m0 <- function(md,p)
{
    md$C <- md.candidates.matrix(md)
    freq <- md %>% group_by(C) %>% count

    function(theta)
    {
        info.mat <- matrix(nrow=p,ncol=p)
        for (j in 1:p)
        {
            for (k in 1:p)
            {
                info.mat[j,k] <- 0.0
                f = freq %>% filter(C[[j]] == T, C[[k]] == T)
                for (i in 1:nrow(f))
                {
                    S <- sum(theta[f$C[i,]])
                    info.mat[j,k] <- info.mat[j,k] + f$n[i] / S^2
                }
            }
        }
        info.mat
    }
}

fisher.scoring <- function(theta0,info,score,eps=1e-2)
{
    repeat
    {
        theta1 <- theta0 + inv(info(theta0)) %*% score(theta0)
        if (max(abs(theta1-theta0)) < eps)
            return(theta1)
        theta0 <- theta1
    }
}



#' maximum likelihood estimator of the
#' parameters of a series system with
#' components that have exponentially
#' distributed lifetimes given a
#' sample of masked data.
#'
#' @param md masked data
#'
#' @return mle of theta given masked data
#' @export
mle.exp.series.m0 = function(md,theta0=NULL)
{
    p <- md.nnodes(md)
    if (is.null(theta0))
        theta0 <- rep(1.0,p)

    fisher.scoring(theta0,
                   info.series.m0(md,p),
                   score.series.m0(md,p))
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
exp.series.observed.info.m0 = function(md)
{
    l = loglike.exp.series.m0(md)
    function(theta)
    {
        hessian(func=l,x=theta)
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
#mle.info.exp.series.m0 = function(n,theta)
#{
#    1
#}

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
#series.exp.mle.cov = function(n,theta)
#{
#    series.exp.mle.info(n,theta)^(-1)
#}

#' bootstrap of the mle's sampling distribution
#'
#' @param masked.data sample of masked data
#' @param r replicates
#'
#' @return mle and bootstrap estimate of its covariance matrix
#' @export
# rseries.exp.mle.cov.bootstrap = function(masked.data,r)
# {
#     n = nrow(masked.data)
#     rate.mle = series.exp.mle(masked.data)
#     rate.bs = NULL
#     for (i in 1:r)
#     {
#         d = rmasked.data(n,rate.mle)
#         rate.bs = rbind(rate.bs,series.exp.mle(d))
#     }
#     list(mle=rate.mle,Sigma=cov(rate.bs))
# }
#
#
#
#
#
#
# mse = function(theta,n) { length(theta)*sum(theta)/n }
#
#
# info <- function(l,n)
# {
#     l1 = l[1]
#     l2 = l[2]
#     l3 = l[3]
#
#     A = matrix(c(1/(l1+l2) + 1/(l1+l3), 1/(l1+l2), 1/(l1+l3),
#                  1/(l1+l2), 1/(l1+l2) + 1/(l2+l3), 1/(l2+l3),
#                  1/(l1+l3), 1/(l2+l3), 1/(l1+l3) + 1/(l2+l3)),byrow=T,nrow=3)
#     n/(2.0*(l1+l2+l3))*A
# }
#
#
# loglike <- function(masked_data)
# {
#     m <- masked_data$num_nodes
#     function(theta)
#     {
#         sum(log(apply(
#             masked_data$c,1,function(r) { sum(theta[r]) }))) -
#             sum(masked_data$s) * sum(theta)
#     }
# }
#
# rmasked_data <- function(n,theta,w)
# {
#     m <- length(theta)
#     t <- matrix(stats::rexp(n*m,rate=theta),ncol=m,byrow=T)
#     k <- integer(length=n)
#     s <- numeric(length=n)
#     for (i in 1:n)
#     {
#         k[i] <- which.min(t[i,])
#         s[i] <- t[i,k[i]]
#     }
#
#     c <- matrix(nrow = n, ncol = m, F)
#
#     if (is.null(w)) {
#         w <- rep(m, n)
#     }
#     w <- c(w, rep(w[length(w)], n - length(w)))
#
#     for (i in 1:n)
#     {
#         c[i, k[i]] <- T
#         c[i, sample((1:m)[-k[i]], size = w[i] - 1, replace = F)] <- T
#     }
#     data <- data.frame(s=s,k=k,t=t,w=w)
#     data$c <- c
#     data
# }
#
#
#
#
#
# n=1000
# theta.star=c(1,3,5)
#
#
#
# data=rmasked_data(n,theta.star,2)
# ll=loglike(data)
# theta.mle=optim(c(1,1,1),function(theta) { -ll(theta) })$par
# theta.mle
#
# H.mle=hessian(func=function(theta) { -ll(theta) },x=theta.mle)
# #I.mle=info(theta.mle,n)
#
# cov.mle=inv(H.mle)
# cov.mle
#
# cov.star=inv(info(theta.star,n))
# cov.star
#
# #delta=theta.mle - theta.star
#
# mse.mle = sum(diag(cov.mle))
# mse.a = mse(theta.mle,n)
# mse.star = mse(theta.star,n)
#
# mse.mle
# mse.star
# mse.a
#
#
#
# #### sufficient statistics
#
# d = rmasked_data(10,theta.star,2)
# loglike_ss = function(masked_data)
# {
#     n = nrow(masked_data)
#     m = ncol(masked_data$c)
#     freq = d %>% group_by(c) %>% count
#
#     #t(apply(freq$c,1,function(r) { (1:m)[r]}))
#
# }
#
# loglike_ss(data)
