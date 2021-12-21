#' Construct a new series system
#'
#' @description
#' Constructs a new series system from the specified
#' list of univariate distribution functions, the
#' pdfs and cdfs. The distribution functions can be
#' a mix of any parametric families.
#'
#' @param pdfs Parametric pdf functions for components.
#'             each should accept a failure time and a
#'             parameter matrix theta, which denotes
#'                 fk(t|theta).
#'
#' @param cdfs Parametric cdf functions for components.
#'             each should accept a failure time and a
#'             parameter matrix theta, which denotes
#'                 Fk(t|theta).
#'
#' @return Object representing the series system.
#' @export
new_series_system = function(pdfs,cdfs)
{
  #stopifnot(is_pdf(pdfs))
  #stopifnot(is_cdf(cdfs))

  pdfs[[1]]$param.dim

  structure(list(pdfs=pdfs,cdfs=cdfs),
            param.dim=c(length(pdfs),p),
            class=c("series_system","random_variable"))
}

new_series_system_pdf <- function(series) {
  # stopifnot(is_series_system(series))

  m <- attributes(series)$param.dim[1]
  p <- attributes(series)$param.dim[2]

  structure(
    function(t, theta) {
      # stopifnot(is.matrix(theta))
      # stopifnot(nrow(theta) == m)
      # stopifnot(ncol(theta) == p)
      theta <- as.matrix(theta, nrow = m, ncol = p)
      d <- 0.
      for (j in 1:m)
      {
        dj <- series$pdfs[[j]](t, theta[j, ])
        for (k in (1:m)[-j]) {
          dj <- dj * (1 - series$cdfs[[k]](t, theta[k, ]))
        }
        d <- d + dj
      }
      d
    },
    param.dim = c(m, p),
    class = c("series_system_pdf", "pdf", "function")
  )
}

new_series_system_cdf = function(series)
{
  # stopifnot(is.series.system(series))

  m = attributes(series)$param.dim[1]
  p = attributes(series)$param.dim[2]

  structure(
    function(t,theta)
    {
      #stopifnot(is.matrix(theta))
      #stopifnot(nrow(theta) == m)
      #stopifnot(ncol(theta) == p)
      theta=as.matrix(theta,nrow=m,ncol=p)
      d = 1.
      for (j in 1:m)
        d = d * (1-series$cdfs[[j]](t,theta[j,]))
      return(1-d)
    },
    param.dim=c(m,p),
    class=c("series_system_cdf","cdf","function")
  )
}

