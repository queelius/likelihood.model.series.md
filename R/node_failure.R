#' Constructs a pdf object for the conditional node failure
#' in an exponential series system according to candidate model
#' m0, \code{f(k|c,s)}.
#'
#' \code{f(k|s,c) = h_k(s)/h(s) I(k in c)}. This simplifies to
#' \code{f(k|c) = theta[k] / sum(theta[j],j in c)} for
#' the exponential series system.
#'
#' @param exp_series An \code{exp_series} object to construct the
#'                   node failure pdf object for.
#' @export
make_exp_series_node_failure_m0 <- function(exp_series)
{
    theta <- param(exp_series)
    m <- num_nodes(exp_series)
    function(k,c,s)
    {
        if (!(k %in% (1:m)[c]))
            0
        else
            theta[k] / sum(theta[c])
    }
}

#' Constructs a pdf object for the conditional node failure
#' in an exponential series system according to candidate model
#' m1, \code{f(k|c,s,alpha)}.
#'
#' \code{f(k|s,c,alpha) = alpha h_k(s)/h(s) I(k in c) + (1-alpha) h_k(s)/h(s) I(k not in c)}.
#' This simplifies to \code{f(k|c,alpha) =
#'     alpha*theta[k]/sum(theta[j],j in c) I(k in c) +
#'     (1-alpha)*theta[k]/sum(theta[j],j in c) I(k not in c)}
#' for the exponential series system.

#' @export
make_exp_series_node_failure_m1 <- function(exp_series)
{
    theta <- param(exp_series)
    m <- num_nodes(exp_series)
    function(k,c,s,alpha)
    {
        if ((k %in% (1:m)[c]))
            theta[k] / sum(theta[c])
    }
}

#' Decorate masked data (tbl_md) with node failure probabilities.
#'
#' Under model m0, we do not know which node caused the failure,
#' assuming \code{|C| > 1}, but if we have an estimate of the node
#' hazard functions (e.g., an estimate of theta from masked data),
#' then we can estimate the probability of node failure.
#'
#' We decorate masked data \code{md} with an estimate of the
#' probabilities, \code{f(k|s,c)} for \code{k=1,...,k=m}.
#'
#' @param md masked data
#' @param h hazard funct
#' @export
md_series_node_failure_decorator_m0 <- function(md,h)
{
    fk <- md_node_failure_m0(h)
    m <- md_num_nodes(md)
    n <- nrow(md)

    # if md does not have a column for candidates,
    # then we assume all nodes are candidates.
    C <- md_candidates_to_matrix(md)

    md$K <- matrix(rep(m*n),ncol=m)
    for (i in 1:n)
    {
        for (k in 1:m)
            md$K[i,k] <- fk(k,md$s[i],C[i,])
    }
    md
}
