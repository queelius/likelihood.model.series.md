#' #' Generic method for obtaining the hazard function of
#' #' an object.
#' #'
#' #' @param x The object to obtain the hazard function of
#' #' @param ... Additional arguments to pass
#' #' @export
#' hazard <- function(x, ...)
#' {
#'     UseMethod("hazard",x)
#' }
#'
#' #' Method for obtaining the hazard function of
#' #' a \code{dist} object.
#' #'
#' #' @param x The object to obtain the hazard function of
#' #' @param ... Additional arguments to pass
#' #' @export
#' hazard.dist <- function(x, ...)
#' {
#'     pdfx <- pdf(x,...)
#'     survx <- surv(x,...)
#'     function(t) pdfx(t)/survx(t)
#' }
#'
#'
#' #' Generic method for obtaining the pdf function of
#' #' an object.
#' #'
#' #' @param x The object to obtain the hazard function of
#' #' @param ... Additional arguments to pass
#' #' @export
#' pdf <- function(x, ...)
#' {
#'     UseMethod("pdf",x)
#' }
#'
#' #' Generic method for obtaining the survival function of
#' #' an object.
#' #'
#' #' @param x The object to obtain the survival function of
#' #' @param ... Additional arguments to pass
#' #' @export
#' surv <- function(x, ...)
#' {
#'     UseMethod("surv",x)
#' }
#'
#' #' Generic method for obtaining the survival function of
#' #' a \code{dist} (distribution) object.
#' #'
#' #' @param x The object to obtain the cdf of
#' #' @param ... Additional arguments to pass
#' #' @export
#' cdf <- function(x, ...)
#' {
#'     UseMethod("cdf",x)
#' }
#'
#' #' Method for obtaining the cdf of a \code{dist} (distribution) object.
#' #'
#' #' @param x The \code{dist} object to obtain the cdf of
#' #' @param ... Additional arguments to pass
#' #' @export
#' cdf.dist <- function(x,...)
#' {
#'     survx <- surv(x,...)
#'     function(t) 1-survx(t)
#' }
#'
#' #' Mmethod for obtaining the survival function a \code{dist} (distribution) object.
#' #'
#' #' @param x The \code{dist} object to obtain the survival function of
#' #' @param ... Additional arguments to pass
#' #' @export
#' surv.dist <- function(x,...)
#' {
#'     cdfx <- cdf(x,...)
#'     function(t) 1-cdfx(t)
#' }
#'
#' #' Method for obtaining the number of components in a system object.
#' #'
#' #' @param x The object to obtain the number of components of
#' #' @param ... Additional arguments to pass
#' #' @export
#' num_comp <- function(x,...)
#' {
#'     UseMethod("num_comp",x)
#' }
#'
#' #' Method for obtaining the number of components in series system object.
#' #'
#' #' @param x The object to obtain the number of components of
#' #' @param ... Additional arguments to pass
#' #' @export
#' num_comp.series <- function(x,...)
#' {
#'     x$num_comp
#' }
