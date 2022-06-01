#' Component failure probability decorator.
#'
#' We decorate masked data \code{md} with an estimate of the
#' probabilities, \code{Pr{K[i]=j|T[i]=t[i],C[i])}}.
#'
#' \code{md} must have columns \code{x1},...,\code{xm} representing candidate
#' sets as boolean vectors, \code{s} or \code{t} respectively representing
#' right-censored series system lifetime, and if right-censoring, then there
#' must also be \code{delta} representing the right-censoring indicator.
#'
#' @param md masked data
#' @param h a vector of hazard functions for the component lifetimes in the
#'          series system. \code{h[j]} should refer to the
#'          hazard function for component indexed by \code{j}.
#' @return \code{md} decorated with component probability failures.
#' @importFrom tibble tibble
#' @importFrom dplyr %>%
#' @importFrom md.tools md_decode_matrix
#' @importFrom dplyr bind_cols
#' @export
md_series_component_failure_probability_decorator <- function(md,h)
{
    #print(md)
    # m <- length(h)
    # C <- md_decode_matrix(md,"x")
    # right_censoring <- "delta" %in% colnames(md)
    # stopifnot(!right_censoring || "s" %in% colnames(md))
    #
    # t <- ifelse(right_censoring, md$s, md$t)
    # n <- nrow(md)
    #
    # pk <- matrix(0,nrow=n,ncol=m)
    # for (i in 1:n)
    # {
    #     if (!right_censoring || (right_censoring && !md$delta[i]))
    #     {
    #         h.C <- function(t)
    #         {
    #             v <- 0
    #             for (j in (1:m)[C[i,]])
    #                 v <- v + h[[j]](t)
    #             v
    #         }
    #
    #         for (j in (1:m)[C[i,]])
    #             pk[i,j] <- h[[j]](t[i]) / h.C(t[i])
    #     }
    # }
    #
    # pk <- tibble::as_tibble(pk)
    # colnames(pk) <- paste0("pk",1:m)
    # md %>%- dplyr::bind_cols(pk)
}



#' Decorate masked data \code{md} with a some measure of the
#' difference between component failure probabilities.
#'
#' @param md masked data
#' @param pk prefix label of component failure probabilities
#' @param diff difference measure, defaults to range length
#' @importFrom md.tools md_decode_matrix
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @export
md_series_component_failure_decorator_diff <- function(
    md,
    pk="pk",
    diff=function(r) range(r))
{
    n <- nrow(md)
    pk <- md_decode_matrix(md,"pk")
    m <- ncol(pk)

    diffs <- numeric(n)
    for (i in 1:n)
        diffs[i] <- diff(pk[i,])
    md %>% mutate(pk.diffs = diffs)
}
