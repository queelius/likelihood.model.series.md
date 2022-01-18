#' Masked data for a 3-out-of-3 (series system) with exponentially distributed nodes.
#'
#' Masked data containing a sample of 10000 system lifetimes and other attributes
#' where the system is parameterized by theta := (3,3,3) and
#' the candidate model is m0.
#'
#' @format A data frame with 10000 rows and 9 variables:
#' \describe{
#'   \item{s}{Real observable variable, system lifetime}
#'   \item{k}{Integer latent variable, the failed node}
#'   \item{w}{Integer observable variable, number of candidates}
#'   \item{t.1-t.3}{Real latent variable, lifetimes of the 3 nodes}
#'   \item{c.1-c.3}{Boolean observable variable, c.j TRUE indicates nodes j is in candidate set}
#' }
#' @source \url{https://github.com/queelius/masked.data/blob/master/data-raw/data4_generator.R}
"data4"
