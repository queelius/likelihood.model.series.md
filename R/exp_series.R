#' Generates masked data for a series system with exponentially distributed
#' nodes and candidate sets according to candidate_model.
#'
#' @param n Integer. The sample size (each row is an observation).
#' @param theta Numeric vector. The jth component has a failure rate theta_j.
#' @param w Integer vector. For the ith observation, generate w_j candidates.
#' @param candidate_model Function that accepts masked data as an argument.
#'                        The candidate model, defaults to md_candidate_m0.
#'                        If set to NULL, then do not generate a candidate set.
#'                        md_mle_exp_series will treat such masked data as
#'                        a sample that includes every node as candidates.
#' @param metadata Boolean. If TRUE writes meta-data for series system to
#'                 attributes of masked data.
#' @return masked data, a data frame of n observations, (s,k,t1,...,tm,c1,...,cm)
#'         where k, t, and c are covariates (or predictors) of s,k,t1,...,tm.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' md_exp_series(n=10,theta=c(1,2,3),w=rep(2,10))
md_exp_series <- function(n,theta,w,candidate_model=md_candidate_m0,metadata=T)
{
    m <- length(theta)

    t <- tibble::as_tibble(matrix(stats::rexp(n*m,rate=theta), ncol=m, byrow=T))
    names(t) <- paste("t",1:m,sep=".")

    md <- tibble::tibble(
        s = apply(t,1,min),
        k = apply(t,1,which.min),
        w = w)
    md <- dplyr::bind_cols(md,t)
    if (!is.null(candidate_model))
        md <- candidate_model(md,m)

    if (metadata)
    {
        args <- md_func_args()
        nodes <- list()
        for (j in 1:m)
            nodes[[j]] <- list(
                "family" = "exponential",
                "index" = theta[j])
        attr(md,"sim") <- list(
            "family" = toString(args[1]),
            "index" = theta,
            "nodes" = nodes,
            "candidate_model" = toString(args["candidate_model"]))
        attr(md,"masked") <- c("k",paste("t",1:m,sep="."))
        attr(md,"m") <- m
        attr(md,"n") <- n
    }
    md
}

#' Kernel log-likelihood for masked data m0 for exponential series system.
#'
#' The log of the kernel of the likelihood function for masked data
#' for a series system with exponentially distributed lifetimes
#' and candidate sets that model m0.
#'
#' This is an optimization of \link[masked.data]{md_kloglike_exp_series_m0_ref}.
#' @export
md_kloglike_exp_series_m0 <- function(md)
{
    #stopifnot(md_is_masked_data(md))
    m <- md_num_nodes(md)
    C <- md_candidates_as_matrix(md)

    function(rate)
    {
        #stopifnot(length(rate) == m)
        sum(log(apply(
            C,1,function(r) { sum(rate[r]) }))) -
            sum(md$s)*sum(rate)
    }
}

#' Kernel log-likelihood for masked data m0 for exponential series system.
#'
#' The log of the kernel of the likelihood function for masked data
#' for a series system with exponentially distributed lifetimes
#' and candidate sets that model m0.
#'
#' This is the unoptimized version, which serves as a ground-truth
#' for testing a more efficient implementation.
#'
#' @param md masked data for candidate model m0
#' @export
md_kloglike_exp_series_m0_ref <- function(md)
{
    C <- md_candidates_as_matrix(md)
    function(rate)
    {
        s <- 0.0
        for (i in 1:nrow(C))
        {
            s <- s + log(sum(rate[C[i,]]))
        }
        s - sum(md$s) * sum(rate)
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
    sum.s = -sum(md$s)
    counts <- md %>% dplyr::group_by_at(dplyr::vars(dplyr::starts_with("c."))) %>% dplyr::count()
    m = ncol(counts)-1

    function(rate)
    {
        v <- matrix(rep(sum.s,m),nrow=m)
        for (j in 1:p)
        {
            f = freq %>% dplyr::filter(paste("c.",j,sep="") == T)
            print(f)
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
    freq <- md %>% dplyr::group_by(C) %>% dplyr::count()
    print(freq)
    p = ncol(md$C)

    function(rate)
    {
        info.mat <- matrix(nrow=p,ncol=p)
        for (j in 1:p)
        {
            for (k in 1:p)
            {
                info.mat[j,k] <- 0.0
                f = freq %>% dplyr::filter(C[[j]]==T,C[[k]]==T)
                print(f)
                for (i in 1:nrow(f))
                {
                    print(i)
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


