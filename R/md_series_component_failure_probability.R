#' Component failure probability decorator.
#'
#' We decorate masked data `md` with an estimate of the
#' probabilities, `Pr{K[i]=j|T[i]=t[i],C[i])}`.
#'
#' `md` must have columns `x1,...,xm` representing candidate
#' sets as boolean vectors, `s` or `t` respectively representing
#' right-censored series system lifetime, and if right-censoring, then there
#' must also be `delta` representing the right-censoring indicator.
#'
#' @param md masked data
#' @param h a vector of hazard functions for the component lifetimes in the
#'          series system. `h[j]` should refer to the
#'          hazard function for component indexed by `j`.
#' @return `md` decorated with component probability failures.
#' @importFrom tibble tibble
#' @importFrom dplyr %>%
#' @importFrom md.tools md_decode_matrix
#' @importFrom dplyr bind_cols
#' @export
md_series_component_failure_probability_C1_C2 <- function(md, h)
{
    m <- length(h)
    C <- md_decode_matrix(md,"x")
    right_censoring <- "delta" %in% colnames(md)
    stopifnot(!right_censoring || "s" %in% colnames(md))

    t <- ifelse(right_censoring, md$s, md$t)
    n <- nrow(md)

    pk <- matrix(rep(0,n*m),nrow=n,ncol=m)
    for (i in 1:n)
    {
        if (!right_censoring || (right_censoring && !md$delta[i]))
        {
            h.C <- function(t) {
                v <- 0
                for (j in (1:m)[C[i,]])
                    v <- v + h[[j]](t)
                v
            }

            for (j in (1:m)[C[i,]])
                pk[i,j] <- h[[j]](t[i]) / h.C(t[i])
        }
    }

    md %>% dplyr::bind_cols(md_encode_matrix(pk, "pk"))
}
