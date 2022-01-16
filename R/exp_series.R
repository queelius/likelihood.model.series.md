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
md_exp_series_m0 <- function(n,rate,w)
{
    m <- length(rate)

    t <- as_tibble(matrix(stats::rexp(n*m,rate=rate), ncol=m, byrow=T))
    names(t) <- paste("t",1:m,sep=".")

    md <- tibble(s = apply(t,1,min),
                 k = apply(t,1,which.min),
                 w = w)
    md <- md_rcandidates_m0(md)
    md <- bind_cols(md,t)

    # set up attribute properties
    nodes <- list()
    for (j in 1:m)
        nodes[[j]] <- list("family" = "exponential",
                           "index"  = rate[j])

    attr(md,"sim") <- list("theta" = rate,
                           "nodes" = nodes,
                           "model" = "m0")
    attr(md,"masked") <- c("k",paste("t",1:m,sep="."))

    md
}

#' @brief log-likelihood function of masked data for a series system
#'        with exponentially distributed lifetimes.
#'
#' @param md masked data for candidate model m0
#'
#' @return log-likelihood function of type R^m -> R
#' @export
md_loglike_exp_series_m0 <- function(md)
{
    m <- md_num_nodes(md)
    C <- md_candidates_as_matrix(md)

    function(rate)
    {
        rate <- as.matrix(rate,nrow=m,ncol=1)
        sum(log(apply(
            C,1,function(r) { sum(rate[r]) }))) -
            sum(md$s) * sum(rate)
    }
}

md_score_exp_series_m0 <- function(md)
{
    md$C <- md_candidates_as_matrix(md)
    freq = md %>% group_by(C) %>% count
    sum.t = -sum(md$s)
    p = ncol(md$C)

    function(rate)
    {
        v <- matrix(rep(sum.t,p),nrow=p)
        for (j in 1:p)
        {
            f = freq %>% filter(C[[j]] == T)
            for (i in 1:nrow(f))
            {
                S <- sum(rate[f$C[i,]])
                v[j,1] <- v[j,1] + f$n[i] / S
            }
        }
        v
    }
}

md_info_exp_series_m0 <- function(md)
{
    md$C <- md_candidates_as_matrix(md)
    freq <- md %>% group_by(C) %>% count
    p = ncol(md$C)

    function(rate)
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
                    S <- sum(rate[f$C[i,]])
                    info.mat[j,k] <- info.mat[j,k] + f$n[i] / S^2
                }
            }
        }
        info.mat
    }
}

md_fisher_scoring <- function(theta0,info,score,eps=1e-2)
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
#' @param rate0 initial guess for MLE
#'
#' @return mle of theta given masked data
#' @export
md_mle_exp_series_m0 = function(md,rate0=NULL)
{
    if (is.null(rate0))
        rate0 <- rep(1.,md_num_nodes(md))

    md_fisher_scoring(rate0,
                      md_info_series_m0(md),
                      md_score_series_m0(md))
}

#' information matrix of the rate parameter
#' of the series system with exponentially
#' distributed component lifetimes given
#' masked data m0.
#'
#' @param md masked data
#'
#' @return observed info
#' @export
md_observed_info_exp_series_m0 = function(md)
{
    l = md_loglike_exp_series_m0(md)
    function(rate)
    {
        hessian(func=l,x=rate)
    }
}


