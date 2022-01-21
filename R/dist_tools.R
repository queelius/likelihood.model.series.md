#' Generic method for obtaining the parameters of
#' a parametric distribution.
#'
#' @param x The object to obtain the parameters of
#'
#' @export
params <- function(x, ...)
{
    UseMethod("params",x)
}

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
num_nodes <- function(series)
{
    series$num_nodes
}

#' Generic method for sampling from distribution objects.
#'
#' @param x The object to sample from.
#'
#' @export
sampler <- function(x, ...)
{
    UseMethod("sampler",x)
}
