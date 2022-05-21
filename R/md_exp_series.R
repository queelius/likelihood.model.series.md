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
md_maked_data_exp_series_dist_reg_cand <- function(n,theta,metadata=T)
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

#' Log-likelihood for masked data with candidate sets obeying the
#' regular candidate set model for exponential series system.
#'
#' @param md masked data
#' @importFrom dplyr %>%
#' @export
md_loglike_exp_series_reg_cand <- function(md)
{
    md$C <- md_decode_matrix(md,"x")
    freq = md %>% dplyr::group_by(C) %>% dplyr::count()
    sum.t = -sum(md$t)
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
#' @param md masked data for regular candidate model
#'
#' @return score function of type R^m -> R
#' @importFrom dplyr %>%
#' @export
md_score_exp_series_reg_cand <- function(md)
{
    t <- -sum(md$t)
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
#' the regular candidate model.
#'
#' @param md masked data with candidate sets \code{x1,...,xm} that meet the
#'           regular candidate model
#' @return observed information matrix of type R^m -> R^(m x m)
#' @importFrom dplyr %>%
#' @export
md_info_exp_series_reg_cand <- function(md)
{
    s <- -sum(md$s)
    md$C <- md_decode_matrix(md,"x")
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


#' Constructs a pdf (pmf) object for the conditional probability of component
#' failure \code{Pr{K[i]=j|C[i]=c[i],T[i]=t[i]) = h_j(t[i])/h(t[i]) I(j in c[i])}.
#' in an exponential series system for the regular candidate model.
#'
#' @param theta parameter value of \code{exp_series_dist}
#' @export
md_exp_series_component_failure_reg_cand <- function(theta)
{
    if (is.matrix(theta))
        theta <- as.vector(theta)

    theta <- unlist(theta)
    m <- length(theta)
    function(k,c,t)
        ifelse(t <= 0 || !(k %in% (1:m)[c]), 0, theta[k] / sum(theta[c]))
}





#' Construct exponential series object.
#'
#' @param rate failure rates
#'
#' @export
make_exp_series_dist <- function(rate)
{
    structure(list(
        theta=unlist(rate),
        num_nodes=length(rate)),
        class=c("exp_series","series_dist","dist"))
}


#' Method for obtaining the variance-covariance of a \code{exp_series_dist} object.
#'
#' @param object The \code{exp_series}The object to obtain the variance of
#' @importFrom stats vcov
#' @export
vcov.exp_series_dist <- function(object, ...)
{
    diag(1/object$theta)^2
}



#' Method to obtain the hazard function of
#' an \code{exp_series_dist} object.
#'
#' @param x The \code{exp_series_dist} object to obtain the hazard function of
#'
#' @export
hazard.exp_series_dist <- function(x, ...)
{
    theta <- params(x)
    function(t,...)
        ifelse(t <= 0,0,sum(theta))
}


#' Method to obtain the pdf of an \code{exp_series_dist} object.
#'
#' @param x The object to obtain the pdf of
#'
#' @export
pdf.exp_series_dist <- function(x, ...)
{
    theta <- params(x)
    function(t,...)
    {
        ifelse(t <= 0,0,sum(theta))
    }
}

#' Method to sample from an \code{exp_series_dist} object.
#'
#' @param x The \code{exp_series_dist} object to sample from.
#' @importFrom algebraic.mle sampler
#' @export
sampler.exp_series_dist <- function(x,...)
{
    rates <- params(x)
    m <- length(rates)
    function(n=1,...)
        series_data(matrix(stats::rexp(m*n,rates,...),nrow=n,ncol=m))
}

#' Method for obtaining the parameters of
#' a \code{series} distribution object.
#'
#' @param x The \code{series} object to obtain the parameters of
#' @importFrom algebraic.mle params
#' @export
params.series <- function(x, ...)
{
    x$theta
}
