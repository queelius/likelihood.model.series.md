#' Masked data for a series system with exponentially distributed nodes
#'
#' A masked dataset containing the system lifetime and other attributes of 100000
#' series system masked data observations with candidate set model m0.
#'
#' @format A data frame with 1000 rows and 9 variables:
#' \describe{
#'   \item{s}{Real observable variable, system lifetime}
#'   \item{k}{Integer latent variable, the failed node}
#'   \item{w}{Integer observable variable, number of candidates}
#'   \item{t.1}{Real latent variable, lifetime of node 1}
#'   \item{t.2}{Real latent variable, lifetime of node 2}
#'   \item{t.3}{Real latent variable, lifetime of node 3}
#'   \item{t.4}{Real latent variable, lifetime of node 4}
#'   \item{t.5}{Real latent variable, lifetime of node 5}
#'   \item{c.1}{Boolean observable variable, TRUE indicates node 1 is in candidate set}
#'   \item{c.2}{Boolean observable variable, TRUE indicates node 2 is in candidate set}
#'   \item{c.3}{Boolean observable variable, TRUE indicates node 3 is in candidate set}
#'   \item{c.4}{Boolean observable variable, TRUE indicates node 1 is in candidate set}
#'   \item{c.5}{Boolean observable variable, TRUE indicates node 2 is in candidate set}
#'
#'   ...
#' }
#' @source \url{https://github.com/queelius/masked.data/blob/master/data-raw/data2_generator.R}
"data2"
