#' A pdf generator of conditional node failure in a series
#' system acording to candidate mdoel m0, \code{f(k|s,c)}.
#'
#' \code{f(k|s,c)} is proportional to \code{h_k(s)/h(s) I(k in c)}.
#'
#' @param h vector of hazard functions for the nodes in series system.
#' @export
md_series_node_failure_m0 <- function(h)
{
    m <- length(h)
    function(k,s,c)
    {
        cs <- (1:m)[c]

        if (!(k %in% cs) || s <= 0)
            return(0.)

        denom <- 0.
        for (j in cs)
            denom <- denom + h[[j]](s)
        h[[k]](s) / denom
    }
}

#' A pdf generator of conditional node failure in a series
#' system acording to candidate mdoel m0, \code{f(k|s,c)}.
#'
#' \code{f(k|s,c,alpha)} is proportional to \code{alpha h_k(s)/h(s) I(k in c) + (1-alpha) h_k(s)/h(s) I(k not in c)}.
#' @param h vector of hazard functions for the nodes in series system.
#' @export
md_series_node_failure_m1 <- function(h)
{
    function(k,s,c,alpha)
    {
        0.
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
