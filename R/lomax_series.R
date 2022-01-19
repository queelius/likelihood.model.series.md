#' Generates masked data for a series system with lomax distributed
#' nodes and candidate sets according to candidate_model.
#'
#' @param n Integer. The sample size (each row is an observation).
#' @param lambda Numeric vector.
#' @param kappa Numeric vector. The jth node is parameterized by theta_j := (lambda_j,kappa_j).
#' @param w Integer vector. For the ith observation, generate w_j candidates.
#' @param candidate_model Function that accepts masked data as an argument.
#'                        The candidate model, defaults to md_candidate_m0.
#'                        If set to NULL, then do not generate a candidate set.
#'                        md_mle_exp_series will treat such masked data as
#'                        a sample that includes every node as candidates.
#' @param metadata Boolean. If TRUE writes meta-data for series system to
#'                 attributes of masked data (tbl_md).
#' @return masked data, a data frame of n observations, (s,k,t1,...,tm,c1,...,cm)
#'         where k, t, and c are covariates (or predictors) of s,k,t1,...,tm.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' md_lomax_series(n=10,lambda=c(1,2,3),kappa=c(4,5,6),w=rep(2,10))
md_lomax_series = function(n,lambda,kappa,w,candidate_model=md_candidate_m0,metadata=T)
{
    m = length(lambda)
    t <- matrix(extraDistr::rlomax(n*m,lambda=lambda,kappa=kappa),ncol=m,byrow=T)
    md <- md_series_data(t,w,candidate_model)

    if (metadata)
    {
        args <- md_func_args()
        nodes <- list()
        for (j in 1:m)
            nodes[[j]] <- list(
                "family" = "lomax",
                "index"  = c(lambda[j],kappa[j]))
        attr(md,"sim") <- list(
            "family" = toString(args[1]),
            "index" = list(lambda=lambda,kappa=kappa),
            "nodes" = nodes,
            "candidate_model" = toString(args["candidate_model"]))
        attr(md,"masked") <- c("k",paste("t",1:m,sep="."))
        attr(md,"m") <- m
        attr(md,"n") <- n
    }
    md
}

#' Kernel log-likelihood for masked data m0 for lomax series system.
#'
#' The log of the kernel of the likelihood function for masked data
#' for a series system with lomax distributed lifetimes
#' and candidate sets that model m0.
#'
#' This is the unoptimized version, which serves as a ground-truth
#' for testing a more efficient implementation.
#'
#' @param md masked data for candidate model m0
#' @export
md_kloglike_lomax_series_m0_ref <- function(md)
{
    C <- md_candidates_as_matrix(md)
    m <- md_num_nodes(md)
    function(theta)
    {
        lambda <- theta[1:m]
        kappa <- theta[(m+1):length(theta)]

        s <- 0.0
        for (i in 1:nrow(C))
        {
            tmp <- 0
            for (k in (1:m)[C[i,]])
            {
                tmp <- tmp + log(lambda[k] * kappa[k] / (1 + lambda[k]))
            }
            s <- s + log(tmp)
            for (j in 1:m)
            {
                s <- s - kappa[j] * log(1+lambda[j]*md$s[i])
            }
        }
        s
    }
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
md_info_lomax_series_m0 = function(md)
{
    l = md_loglike_lomax_series_m0(md)
    function(theta)
    {
        -numDeriv::hessian(func=l,x=theta)
    }
}
