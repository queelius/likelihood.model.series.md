#' Generates masked data for a series system with exponentially distributed
#' nodes and candidate sets according to candidate_model.
#'
#' @param n Integer. The sample size (each row is an observation).
#' @param theta Numeric vector. The jth component has a failure rate \code{theta[j]}.
#' @param w Integer vector. For the ith observation, generate \code{w[j]} candidates.
#' @param candidate_model Function that accepts masked data as an argument.
#'                        The candidate model, defaults to \code{md_candidate_m0}.
#'                        If set to \code{NULL}, then do not generate a candidate set.
#'                        \code{md_mle_exp_series} will treat such masked data as
#'                        a sample that includes every node as candidates.
#' @param metadata Boolean. If \code{TRUE} writes meta-data for series system to
#'                 attributes of masked data.
#' @return masked data, a data frame of n observations, \code{(s,k,t1,...,tm,c1,...,cm)}
#'         where \code{k}, \code{t}, and \code{c} are covariates (or predictors) of
#'         \code{s,k,t1,...,tm}.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' md_exp_series(n=10,theta=c(1,2,3),w=rep(2,10))
md_exp_series <- function(n,theta,w,candidate_model=md_candidate_m0,metadata=T)
{
    m <- length(theta)
    t <- matrix(stats::rexp(n*m,rate=theta), ncol=m, byrow=T)
    md <- md_series_data(t,w,candidate_model)

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

#' Kernel log-likelihood for masked data m0 for exponential series system
#' using sufficient statistics.
#'
#' The log of the kernel of the likelihood function for masked data
#' for a series system with exponentially distributed lifetimes
#' and candidate sets that model m0 using sufficient statistics.
#'
#' @param md masked data
#' @importFrom dplyr %>%
#' @export
md_kloglike_exp_series_m0 <- function(md)
{
    md$C <- md_candidates_as_matrix(md)
    freq = md %>% dplyr::group_by(C) %>% dplyr::count()
    sum.t = -sum(md$s)
    n = nrow(freq)
    function(theta)
    {
        f <- sum.t * sum(theta)
        for (i in 1:n)
            f <- f + freq$n[i] * log(sum(theta[freq$C[i,]]))
        f
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
    s <- -sum(md$s)
    md$C <- md_candidates_as_matrix(md)
    cnt <- md %>% dplyr::group_by(C) %>% dplyr::count()
    m <- ncol(md$C)

    function(rate)
    {
        v <- rep(s,m)
        for (j in 1:m)
        {
            for (i in 1:nrow(cnt))
            {
                if (cnt$C[i,j])
                    v[j] <- v[j] + cnt$n[i] / sum(rate[cnt$C[i,]])
            }
        }
        v
    }
}

#' Information matrix (observed) for rate parameter
#' with respect to masked data of a series system
#' with exponentially distributed lifetimes and
#' candidate model m0.
#'
#' @param md masked data for candidate model m0
#'
#' @return observed information matrix of type R^m -> R^(m x m)
#' @importFrom dplyr %>%
#' @export
md_info_exp_series_m0 <- function(md)
{
    s <- -sum(md$s)
    md$C <- md_candidates_as_matrix(md)
    cnt <- md %>% dplyr::group_by(C) %>% dplyr::count()
    m <- ncol(md$C)

    function(rate)
    {
        nfo <- matrix(rep(0.0,m*m),nrow=m)
        for (j in 1:m)
        {
            for (k in 1:m)
            {
                for (i in 1:nrow(cnt))
                {
                    if (cnt$C[i,j] && cnt$C[i,k])
                        nfo[j,k] <- nfo[j,k] + cnt$n[i] / sum(rate[cnt$C[i,]])^2
                }
            }
        }
        nfo
    }
}

#' Maximum likelihood estimator of the parameters of a series
#' system with nodes that have exponentially distributed
#' lifetimes given a sample of masked data according to
#' candidate model m0.
#'
#' @param md masked data
#' @param theta0 initial guess for MLE
#' @param eps stopping condition
#' @param max_iterations stop if iterations reaches max_iterations.
#'
#' @return MLE estimate
#' @export
md_mle_exp_series_m0 = function(md,theta0=NULL,eps=1e-5,max_iterations=250L)
{
    if (is.null(theta0))
        theta0 <- rep(1.,md_num_nodes(md))

    res <- md_fisher_scoring(
        theta0,
        md_info_exp_series_m0(md),
        md_score_exp_series_m0(md),
        eps,
        max_iterations)

    structure(list(
        theta.hat=res$theta.hat,
        iterations=res$iterations,
        max_iterations=max_iterations,
        eps=eps,
        score=res$s,
        info=res$info,
        sigma=res$sigma),
        class=c("md_estimate"),
        attributes=list("candidate_model" = "m0"))
}

#' Constructs a pdf object for the conditional node failure
#' in an exponential series system according to candidate model
#' m0, \code{f(k|c,s) = h_k(s)/h(s) I(k in c)}.
#'
#' This simplifies to \code{f(k|c) = theta[k] / sum(theta[j],j in c)} for
#' the exponential series system.
#'
#' @param theta parameter value of \code{exp_series}
#' @export
md_exp_series_node_failure_m0 <- function(theta)
{
    if (is.matrix(theta))
        theta <- as.vector(theta)

    theta <- unlist(theta)
    m <- length(theta)
    function(k,c,s)
    {
        if (s <= 0 || !(k %in% (1:m)[c]))
            0
        else
            theta[k] / sum(theta[c])
    }
}


#' Constructs the shortest interval for the system lifetime
#' given a candidate set under model m0 with a probability \code{p}
#' that the interval contains the system failure.
#'
#' @param theta parameter value of \code{exp_series}
#' @param p probability that system failure time is in the computed interval
#' @export
md_exp_series_system_failure_interval_m0 <- function(theta,p)
{
    if (is.matrix(theta))
        theta <- as.vector(theta)
    theta <- unlist(theta)
    function(c)
    {
        c(0,qexp(p,sum(theta[c])))
    }
}
