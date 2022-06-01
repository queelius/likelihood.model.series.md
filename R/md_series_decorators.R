#' Generates failure times and component cause of failure for a series
#' system with the given data frame of component lifetimes encoded by
#' the matrix columns prefixed with \code{t} in masked data frame \code{md}.
#'
#' @param md a data frame with the indicated component times-to-failure
#' @return masked data decorated with series lifetime and component cause of failure
#' @importFrom md.tools md_decode_matrix
#' @importFrom md.tools md_mark_latent
#' @importFrom dplyr %>%
#' @export
md_series_lifetime <- function(md)
{
    t <- md_decode_matrix(md,"t")
    md$k <- apply(t,1,which.min)
    md$t <- apply(t,1,min)
    md %>% md_mark_latent("k")
}

#' Generates right-censored system failure times \code{s} and right-censoring
#' indicators \code{delta} for a series system with the given data frame of
#' component lifetimes encoded by the matrix columns prefixed with \code{t}
#'
#' @param md a data frame with the indicated component times-to-failure
#' @param tau vector of right-censoring times
#' @return masked data decorated with series lifetime and component cause of failure
#'
#' @importFrom md.tools md_decode_matrix
#' @importFrom md.tools md_mark_latent
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @export

md_series_lifetime_right_censoring <- function(md,tau)
{
    md %>% mutate(tau = tau) %>%
        mutate(s = ifelse(tau < t, tau, t)) %>%
        mutate(delta = t > tau) %>%
        md_mark_latent("t")
}
