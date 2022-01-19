#' Masked data for a 3-out-of-3 (series system) with exponentially distributed nodes.
#'
#' Masked data containing a sample of 10000 system lifetimes and other attributes
#' where the system is parameterized by \code{theta=c(1,1,1)} and
#' candidate model m0.
#'
#' Candidate set sizes are \code{w=2}.
#'
#' @format A data frame with 10000 rows and 9 variables:
#' \describe{
#'   \item{s}{Real observable variable, system lifetime}
#'   \item{k}{Integer latent variable, the failed node}
#'   \item{w}{Integer observable variable, number of candidates}
#'   \item{t.1}{Real latent variable, lifetime of node 1}
#'   \item{t.2}{Real latent variable, lifetime of node 2}
#'   \item{t.3}{Real latent variable, lifetime of node 3}
#'   \item{c.1}{Boolean observable variable, c.1 TRUE indicates nodes j is in candidate set}
#'   \item{c.2}{Boolean observable variable, c.2 TRUE indicates nodes j is in candidate set}
#'   \item{c.3}{Boolean observable variable, c.3 TRUE indicates nodes j is in candidate set}
#' }
#' @source \url{https://github.com/queelius/masked.data/blob/master/data-raw/exp_series_data_4_gen.R}
"exp_series_data_4"
