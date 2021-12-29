#' Construct a new distribution
#'
#' @param pdf pdf function.
#' @param cdf cdf function.
#'
#' @return Distribution object
#' @export
new_dist <- function(pdf, cdf) {
  stopifnot(is_pdf(pdf))
  stopifnot(is_cdf(cdf))
  stopifnot(all(domain_dim(pdf) == domain_dim(cdf)))
  if (is_parametric(pdf) | is_parametric(cdf)) {
    stopifnot(is_parametric(pdf) == is_parametric(cdf))
    stopifnot(all(param_dim(pdf) == param_dim(cdf)))
  }

  structure(
    list(
      pdf = pdf,
      cdf = cdf
    ),
    domain.dim = domain_dim(pdf),
    codomain.dim = codomain_dim(pdf),
    param.dim = param_dim(pdf),
    class = c("distribution")
  )
}

domain_dim <- function(f) {
  f$domain.dim
}

codomain_dim <- function(f) {
  f$codomain.dim
}

param_dim <- function(f) {
  f$param.dim
}

is_pdf <- function(f) {
  inherits(f, "pdf")
}

is_cdf <- function(f) {
  inherits(f, "cdf")
}

is_parametric <- function(f) {
  !is.null(f$param.dim)
}

domain_dim <- function(f) {
    ifelse(is_parametric(f),length(formals(f))-1, length(formals(f)))
}

new_parametric_pdf <- function(pdf, r, c) {
  stopifnot(is.function(f))

  structure(
    function(t, theta) {
      theta <- as.matrix(theta, nrow = r, ncol = c)
      pdf(t, theta)
    },
    domain.dim = c(1, 1),
    codomain.dim = c(1, 1),
    param.dim = c(r, c),
    class = c("pdf", "function")
  )
}

new_parametric_cdf <- function(cdf, r, c) {
  stopifnot(is.function(f))

  structure(
    function(t, theta) {
      theta <- as.matrix(theta, nrow = r, ncol = c)
      cdf(t, theta)
    },
    domain.dim = c(1, 1),
    codomain.dim = c(1, 1),
    param.dim = c(r, c),
    class = c("cdf", "function")
  )
}
