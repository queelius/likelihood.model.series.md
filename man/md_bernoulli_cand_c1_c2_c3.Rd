% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/md_series_masked_component_cause.R
\name{md_bernoulli_cand_c1_c2_c3}
\alias{md_bernoulli_cand_c1_c2_c3}
\title{Bernoulli candidate set model for systems with unobserved components.}
\usage{
md_bernoulli_cand_c1_c2_c3(
  df,
  p,
  prob = "q",
  comp = "t",
  right_censoring_indicator = "delta"
)
}
\arguments{
\item{df}{masked data.}

\item{p}{a vector of probabilities (p\link{j} is the probability that the jth
system will include a non-failed component in its candidate set,
assuming the jth system is not right-censored).}

\item{prob}{column prefix for component probabilities, defaults to
\code{q}, e.g., \verb{q1, q2, q3}.}

\item{comp}{column prefix for component lifetimes, defaults to \code{t},
e.g., \verb{t1, t2, t3}.}

\item{right_censoring_indicator}{right-censoring indicator column name.
if TRUE, then the system lifetime is right-censored, otherwise it is
observed. If NULL, then no right-censoring is assumed. Defaults to
\code{delta}.}
}
\description{
Bernoulli candidate set model is a particular type of \emph{uninformed} model.
Note that we do not generate candidate sets with this function. See
\code{md_cand_sampler} for that.
}
\details{
This model satisfies conditions C1, C2, and C3.
The failed component will be in the corresponding candidate set with
probability 1, and the remaining components will be in the candidate set
with probability \code{p} (the same probability for each component). \code{p}
may be different for each system, but it is assumed to be the same for
each component within a system, so \code{p} can be a vector such that the
length of \code{p} is the number of systems in the data set (with recycling
if necessary).
}
