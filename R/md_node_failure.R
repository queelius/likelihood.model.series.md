#' Decorate masked data (tbl_md) with node failure probabilities.
#'
#' Under model m0, we do not know which node caused the failure,
#' (note: if |C|=1, under m0 we know precisely which node failed),
#' but if we have an estimate (or know) theta, then we may
#' construct f(k|s,c) and compute the node failure probabilities
#' in a masked data object \code{md}.
#'
#' We decorate masked data \code{md} with an estimate of the
#' probabilities, \code{f(k|s,c)} for \code{k=1,...,k=m}
#' and return the result.
#'
#' @param md masked data
#' @param fk pdf f(k|s,c)
#'
#' @export
md_series_node_failure_decorator_m0 <- function(md,fk)
{
    m <- md_num_nodes(md)
    n <- nrow(md)

    C <- md_candidates_as_matrix(md)
    K <- matrix(rep(NA,m*n),ncol=m)
    for (i in 1:n)
    {
        for (k in 1:m)
            K[i,k] <- fk(k,C[i,],md$s[i])
    }

    K <- as_tibble(K)
    colnames(K) <- paste0("k.",1:m)
    dplyr::bind_cols(md,K)
}
