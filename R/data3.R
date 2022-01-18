#' Masked data for a 10-out-of-10 (series system) with exponentially distributed nodes.
#'
#' Masked data containing a sample of 100000 system lifetimes and other attributes
#' where the system is parameterized by theta := (3,5,4,6,7,2,8,9,10,11) and
#' the candidate model is m0.
#'
#' @format A data frame with 100000 rows and 23 variables:
#' \describe{
#'   \item{s}{Real observable variable, system lifetime}
#'   \item{k}{Integer latent variable, the failed node}
#'   \item{w}{Integer observable variable, number of candidates}
#'   \item{t.1-t.10}{Real latent variable, lifetimes of the 10 nodes}
#'   \item{c.1-c.10}{Boolean observable variable, c.j TRUE indicates nodes j is in candidate set}
#' }
#' @source \url{https://github.com/queelius/masked.data/blob/master/data-raw/data3_generator.R}
"data3"
