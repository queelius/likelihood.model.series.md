% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exp_series_md_c1_c2_c3.R
\name{score.exp_series_md_c1_c2_c3}
\alias{score.exp_series_md_c1_c2_c3}
\title{Generates a score function for an exponential series system (or
competing risks) with respect to parameter \code{theta} for masked component
failure (or masked competing risks) with candidate sets (competing risks)
that satisfy conditions C1, C2, and C3.}
\usage{
score.exp_series_md_c1_c2_c3(model)
}
\arguments{
\item{model}{the exponential series system masked data model (competing risks)
for candidate set model that satisfies C1, C2, and C3}
}
\value{
score function that takes the following arguments:
\itemize{
\item \code{df}: masked data
\item \code{theta}: rate / scale parameters of component lifetime distributions
\item \code{sys.var}: system lifetime (optionally right-censored) column name,
defaults to \code{t}
\item \code{candset}: prefix of Boolean matrix encoding of candidate sets, defaults
to \code{x}, e.g., \verb{x1,...,xm}.
}
}
\description{
Generates a score function for an exponential series system (or
competing risks) with respect to parameter \code{theta} for masked component
failure (or masked competing risks) with candidate sets (competing risks)
that satisfy conditions C1, C2, and C3.
}
