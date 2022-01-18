#' Masked data for a series system with lomax distributed nodes and candidate
#' sets that model m0
#'
#' A masked dataset containing the system lifetime and other attributes of 10000
#' series system with parameters lambda=(3,4,5) and kappa=(2,3,4). Every candidate
#' set (of model m0) has w=2 candidate nodes.
#'
#' @format A data frame with 10000 rows and 9 variables:
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
#'   ...
#' }
#' @source \url{https://github.com/queelius/masked.data/blob/master/data-raw/data1_lomax.R}
"data1_lomax"
