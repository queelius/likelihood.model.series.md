#' Decorate masked data (tbl_md) with component failure probabilities.
#'
#' Under conditions \code{Pr{K[i] in C[i]}} and
#' \code{Pr{C[i] = c[i] | K[i] = j, T[i] = t[i]} = Pr{C[i] = c[i] | K[i] = j', T[i] = t[i]}}
#' for any \code{j in C[i]} and \code{j' in C[i]}, we do not know which component
#' caused the series system failure (unless \code{|C|=1}), but if we have an
#' estimate (or know) \code{theta}, then we may compute \code{Pr{K[i]=j|T[i]=t[i],C[i])}, the component failure
#' probabilities, for a masked data object \code{md}.
#'
#' We decorate masked data \code{md} with an estimate of the
#' probabilities, \code{Pr{K[i]=j|T[i]=t[i],C[i])}},
#' \code{Pr{K[i]=j|T[i]=t[i]}}, and \code{Pr{K[i]=j}} for \code{j in {1,...,m}}.
#'
#' @param md masked data
#' @param Pk probability \code{Pr{K[i]=j|...}} for \code{j in {1,...,m}},
#'
#' @export
md_series_component_failure_decorator <- function(md,Pk,var.t="t",var.cand="x")
{
    m <- md_num_nodes(md)
    n <- nrow(md)

    C <- md_decode_matrix(md,var.cand)
    K <- matrix(rep(NA,m*n),ncol=m)
    for (i in 1:n)
    {
        if (is.na(md[var.t]) && is.na(var.cand))
        {
            for (j in 1:m)
                K[i,j] <- Pk(j)
        }
        else if (is.na(md[var.cand]))
        {
            for (j in 1:m)
                K[i,j] <- Pk(k,C[i,],md[i,var.t])
        }
        else if (is.na(var.t))
        {
            for (j in 1:m)
                K[i,j] <- Pk(k,C[i,],md[i,var.t])
        }

    }

    K <- as_tibble(K)
    colnames(K) <- paste0("pk",1:m)
    dplyr::bind_cols(md,K)
}


#' Decorate masked data \code{md} with a some measure of the
#' difference between component failure probabilities.
#'
#' @param md masked data
#' @param pk prefix label of component failure probabilities
#' @param diff difference measure, defaults to range length
#'
#' @export
md_series_component_failure_decorator_diff <- function(
    md,
    pk="pk",
    diff=function(r) range(r))
{
    m <- md_num_nodes(md)
    n <- nrow(md)

    pk <- md_decode_matrix(md,"k")
    diffs <- numeric(n)
    for (i in 1:n)
        diffs[i] <- q(pk[i,])
    md %>% mutate(pk.diffs = diffs)
}
