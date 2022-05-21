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
#' @param x The object to obtain the number of components of
#'
#' @export
md_num_comp <- function(x)
{
    x$num_comp
}

#' Generates system failure time and component failure for a series system with
#' the given matix of component failure times.
#'
#' @param t matrix of node failure times
#' @importFrom dplyr %>%
#' @export
series_data <- function(t)
{
    m <- ncol(t)
    data <- tibble::tibble(
        k = apply(t,1,which.min),
        t = apply(t,1,min))

    t <- tibble::as_tibble(t)
    names(t) <- paste0("t",1:m)

    data %>% dplyr::bind_cols(t)
}
