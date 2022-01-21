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
#' @param q interval computer for s|c
#'
#' @export
md_series_system_failure_decorator_m0 <- function(md,q)
{
    m <- md_num_nodes(md)
    n <- nrow(md)

    # if md does not have a column for candidates,
    # then we assume all nodes are candidates.
    C <- md_candidates_as_matrix(md)
    lower <- numeric(n)
    upper <- numeric(n)
    for (i in 1:n)
    {
        res <- q(C[i,])
        lower[i] <- res[1]
        upper[i] <- res[2]
    }

    md %>% mutate(s.lower = lower, s.upper = upper)
}
