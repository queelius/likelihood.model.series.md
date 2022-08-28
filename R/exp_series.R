#' Generates a log-likelihood for an exponential series system with respect to
#' parameter \code{theta} for masked data with candidate sets that approximately
#' satisfy conditions C1, C2, and C3.
#'
#' @param md masked data
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr count
#' @importFrom dplyr filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_loglike_exp_series_C1_C2_C3 <- function(md)
{
    sum.t <- NULL
    if ("delta" %in% colnames(md))
    {
        stopifnot("s" %in% colnames(md))
        sum.t <- sum(md$s)
        md <- md %>% dplyr::filter(.data$delta==F)
    }
    else
    {
        stopifnot("t" %in% colnames(md))
        sum.t <- sum(md$t)
    }
    md$C <- md.tools::md_decode_matrix(md,"x")

    md <- md %>% dplyr::group_by(.data$C) %>% dplyr::count()
    n <- nrow(md)

    function(theta)
    {
        f <- -sum.t * sum(theta)
        for (i in 1:n)
            f <- f + md$n[i] * log(sum(theta[md$C[i,]]))
        f
    }
}

#' Generates a score function  for an exponential series system with respect to
#' parameter \code{theta} for masked data with candidate sets that approximately
#' satisfy conditions C1, C2, and C3.
#'
#' @param md masked data for regular candidate model
#' @return score function of type \code{R^m -> R}
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr count
#' @importFrom dplyr filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_score_exp_series_C1_C2_C3 <- function(md)
{
    sum.t <- NULL
    if ("delta" %in% colnames(md))
    {
        stopifnot("s" %in% colnames(md))
        sum.t <- sum(md$s)
        md <- md %>% dplyr::filter(.data$delta==F)
    }
    else
    {
        stopifnot("t" %in% colnames(md))
        sum.t <- sum(md$t)
    }
    md$C <- md.tools::md_decode_matrix(md,"x")
    md <- md %>% dplyr::group_by(.data$C) %>% dplyr::count()
    m <- ncol(md$C)

    function(theta)
    {
        v <- rep(-sum.t,m)
        for (j in 1:m)
        {
            for (i in 1:nrow(md))
            {
                if (md$C[i,j])
                    v[j] <- v[j] + md$n[i] / sum(theta[md$C[i,]])
            }
        }
        as.matrix(v)
    }
}

#' Generates the observed information matrix for an exponential series system
#' with respect to parameter \code{theta} for masked data with candidate sets
#' that approximately satisfy conditions C1, C2, and C3.
#'
#' @param md masked data with candidate sets \code{x1,...,xm} that meet the
#'           regular candidate model
#' @return observed information matrix of type \code{R^m -> R^(m x m)}
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr count
#' @importFrom dplyr filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_info_exp_series_C1_C2_C3 <- function(md)
{
    md <- md %>% dplyr::filter(.data$delta==F)
    md$C <- md.tools::md_decode_matrix(md,"x")
    md <- md %>% dplyr::group_by(.data$C) %>% dplyr::count()
    m <- ncol(md$C)
    n <- nrow(md)

    function(theta)
    {
        nfo <- matrix(rep(0,m*m),nrow=m)
        for (j in 1:m)
        {
            for (k in 1:m)
            {
                for (i in 1:nrow(md))
                {
                    if (md$C[i,j] && md$C[i,k])
                        nfo[j,k] <- nfo[j,k] + md$n[i] / sum(theta[md$C[i,]])^2
                }
            }
        }
        nfo
    }
}

#' Qunatile function for exponential series.
#'
#' @param p quantiles
#' @param rates rate parameters for exponential component lifetimes
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X<=x]}, otherwise, \code{P[X > x]}.
#' @param log.p return the log of the survival function
#' @importFrom stats qexp
#' @export
qexp_series <- function(p,rates,lower.tail=T,log.p=F)
{
    qexp(p,sum(rates),lower.tail,log.p)
}

#' Survival function for exponential series.
#'
#' @param n sample size
#' @param rates rate parameters for exponential component lifetimes
#' @importFrom stats rexp
#' @export
rexp_series <- function(n,rates)
{
    rexp(n,sum(rates))
}

#' pdf for exponential series.
#'
#' @param t series system lifetime
#' @param rates rate parameters for exponential component lifetimes
#' @param log return the log of the pdf
#' @importFrom stats dexp
#' @export
dexp_series <- function(t,rates,log=F)
{
    dexp(t,sum(rates),log)
}

#' cdf for exponential series.
#'
#' @param t series system lifetime
#' @param rates rate parameters for exponential component lifetimes
#' @param log.p return the log of the cdf
#' @param lower.tail logical; logical; if TRUE (default), probabilities are \code{P[X<=x]}, otherwise, \code{P[X > x]}.
#' @importFrom stats pexp
#' @export
pexp_series <- function(t,rates,lower.tail=T,log.p=F)
{
    pexp(t,sum(rates),lower.tail,log.p)
}

#' Survival function for exponential series.
#'
#' @param t series system lifetime
#' @param rates rate parameters for exponential component lifetimes
#' @param log.p return the log of the survival function
#' @importFrom stats pexp
#' @export
survival_exp_series <- function(t,rates,log.p=F)
{
    pexp(t,sum(rates),TRUE,log.p)
}

#' Hazard function for exponential series.
#'
#' @param t series system lifetime
#' @param rates rate parameters for exponential component lifetimes
#' @param log.p return the log of the hazard function
#' @export
hazard_exp_series <- Vectorize(function(t,rates,log.p=F)
{
    ifelse(log.p,
           ifelse(t < 0, -Inf, sum(rates)),
           ifelse(t < 0, 0,sum(rates)))
},vectorize.args = "t")
