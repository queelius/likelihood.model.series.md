#' Generic method for obtaining the hazard function of
#' a random variable.
#'
#' @param x The object to obtain the hazard function of
#'
#' @export
hazard <- function(x, ...)
{
    UseMethod("hazard",x)
}

#' Generic method for obtaining the pdf function of
#' a random variable.
#'
#' @param x The object to obtain the hazard function of
#'
#' @export
pdf <- function(x, ...)
{
    UseMethod("pdf",x)
}


#' Method for obtaining the number of nodes in an object.
#'
#' @param series The object to obtain the number of nodes of
#'
#' @export
md_num_comp <- function(series)
{
    series$num_comp
}

