#' Masked data for a series system with exponentially distributed components
#'
#' Masked data containing the system lifetime and other attributes of 1000
#' series system with parameter value \code{theta=c(3,4,5)} and candidate
#' model m0.
#'
#' Each candidate is of size \code{w=2}.
#'
#' @format A data frame with 1000 rows and 9 variables:
#' \describe{
#'   \item{s}{Real observable variable, system lifetime}
#'   \item{k}{Integer latent variable, the failed node}
#'   \item{w}{Integer observable variable, number of candidates}
#'   \item{t.1}{Real latent variable, lifetime of node 1}
#'   \item{t.2}{Real latent variable, lifetime of node 2}
#'   \item{t.3}{Real latent variable, lifetime of node 3}
#'   \item{c.1}{Boolean observable variable, TRUE indicates node 1 is in candidate set}
#'   \item{c.2}{Boolean observable variable, TRUE indicates node 2 is in candidate set}
#'   \item{c.3}{Boolean observable variable, TRUE indicates node 3 is in candidate set}
#' }
#' @source \url{https://github.com/queelius/masked.data/blob/master/data-raw/exp_series_data_1_gen.R}
"exp_series_data_1"
