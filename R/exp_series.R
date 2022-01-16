#' Generates masked data for a series system with exponentially distributed
#' nodes and candidate sets according to candidate_model.
#'
#' @param n Integer. The sample size (each row is an observation).
#' @param theta Numeric vector. The j-th component has a failure rate theta_j.
#' @param w Integer vector. For the ith observation, generate w_j candidates.
#' @param candidate_model the candidate model, defaults to md_candidate_m0.
#'
#' @return masked data, a data frame of n observations, (s,k,t1,...,tm,c1,...,cm)
#'         where k, t, and c are covariates (or predictors) of s,k,t1,...,tm.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' md <- md_exp_series(
#'     n=10,
#'     theta=c(1,2,3),
#'     w=rep(2,10),
#'     candidate_model=md_candidate_m1)
md_exp_series <- function(n,theta,w,candidate_model=md_candidate_m0)
{
    m <- length(theta)

    t <- tibble::as_tibble(matrix(stats::rexp(n*m,rate=theta), ncol=m, byrow=T))
    names(t) <- paste("t",1:m,sep=".")

    md <- tibble::tibble(s = apply(t,1,min),
                 k = apply(t,1,which.min),
                 w = w)
    md <- candidate_model(md)
    md <- dplyr::bind_cols(md,t)

    # set up attribute properties
    nodes <- list()
    for (j in 1:m)
        nodes[[j]] <- list(
            "family" = "exponential",
            "index" = theta[j])
    attr(md,"sim") <- list(
        "family" = "exponential_series",
        "index" = theta,
        "nodes" = nodes,
        "candidate_model" = toString(match.call()[5]))
    attr(md,"masked") <- c("k",paste("t",1:m,sep="."))

    md
}

#' log-likelihood function of masked data for a series system
#' with exponentially distributed lifetimes.
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

#' score function of masked data for a series system
#' with exponentially distributed lifetimes.
#'
#' @param md masked data for candidate model m0
#'
#' @return score function of type R^m -> R
#' @importFrom dplyr %>%
#' @export
md_score_exp_series_m0 <- function(md)
{
    md$C <- md_candidates_as_matrix(md)
    freq = md %>% dplyr::group_by(C) %>% dplyr::count
    sum.t = -sum(md$s)
    p = ncol(md$C)

    function(rate)
    {
        v <- matrix(rep(sum.t,p),nrow=p)
        for (j in 1:p)
        {
            f = freq %>% dplyr::filter(C[[j]] == T)
            for (i in 1:nrow(f))
            {
                S <- sum(rate[f$C[i,]])
                v[j,1] <- v[j,1] + f$n[i] / S
            }
        }
        v
    }
}

#' Expected information matrix for rate parameter
#' with respect to masked data of a series system
#' with exponentially distributed lifetimes and
#' candidate model m0.
#'
#' @param md masked data for candidate model m0
#'
#' @return expected information function of type R^m -> R^(m x m)
#' @importFrom dplyr %>%
#' @export
md_info_exp_series_m0 <- function(md)
{
    md$C <- md_candidates_as_matrix(md)
    freq <- md %>% dplyr::group_by(C) %>% dplyr::count
    p = ncol(md$C)

    function(rate)
    {
        info.mat <- matrix(nrow=p,ncol=p)
        for (j in 1:p)
        {
            for (k in 1:p)
            {
                info.mat[j,k] <- 0.0
                f = freq %>% dplyr::filter(md$C[[j]] == T, md$C[[k]] == T)
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

#' maximum likelihood estimator of the
#' parameters of a series system with
#' components that have exponentially
#' distributed lifetimes given a
#' sample of masked data.
#'
#' @param md masked data
#' @param theta0 initial guess for MLE
#' @param eps stopping condition
#'
#' @return mle of theta given the observed masked data md
#' @export
md_mle_exp_series_m0 = function(md,theta0=NULL,eps=1e-5)
{
    if (is.null(theta0))
        theta0 <- rep(1.,md_num_nodes(md))

    md_fisher_scoring(theta0,
                      md_info_exp_series_m0(md),
                      md_score_exp_series_m0(md),
                      eps)
}

#' Observed information matrix of the rate parameter
#' of the series system with exponentially
#' distributed component lifetimes given
#' masked data with candidate sets according to model m0.
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
        numDeriv::hessian(func=l,x=rate)
    }
}


