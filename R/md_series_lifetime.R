#' Generates right-censored system failure times and right-censoring
#' indicators for a series system with the given data frame of
#' component lifetimes.
#'
#' @param md a data frame with the indicated component times-to-failure
#' @param tau vector of right-censoring times, defaults to
#'            `Inf` (no right censoring)
#' @param compvar component lifetimes
#' @param sysvar right-censored system lifetime
#' @param delta right-censoring indicator, defaults to "delta"
#' @return masked data with right-censoring related columns
#'
#' @importFrom md.tools md_decode_matrix md_mark_latent
#' @importFrom dplyr %>%
#' @export
md_series_lifetime_right_censoring <- function(md,
                                               tau=Inf,
                                               compvar="t",
                                               sysvar=NULL,
                                               delta="delta")
{
    # retrieve component lifetimes as a matrix
    ts <- md_decode_matrix(md, compvar)
    stopifnot(!is.null(ts))

    if (is.null(sysvar))
        sysvar <- compvar

    sys <- apply(ts, 1, min)
    md[, sysvar] <- ifelse(tau < sys, tau, sys)
    md[, delta] <- sys > tau
    md %>% md_mark_latent(paste0(compvar, 1:(ncol(ts))))
}
