#' Likelihood model for exponential series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3,
#' desribed below. 
#' 
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood_model` package by providing the following methods:
#' 
#'  (1) `loglik.exp_series_md_c1_c2_c3`
#'  (2) `score.exp_series_md_c1_c2_c3`
#'  (3) `hess_loglik.exp_system_series_md_c1_c2_c3`
#' 
#' These are useful for doing maximum likelihood estimation, hypothesis
#' testing (e.g., likelihood ratio test), estimation of asymptotic sampling
#' distribution given data from the DGP according to the specified model,
#' etc.
#' 
#' It is designed to work well with the `likelihood_model` R package. In
#' particular, it is intended to be used with the `likelihood_contr_model`
#' object, which is a `likelihood_model` object that allows likelihood
#' contributions to be added for whatever data model you have in mind.
#' 
#' In this likelihood model, masked component data approximately satisfies the
#' following conditions:
#' 
#' C1: `Pr{K[i] in C[i]) = 1`
#' C2: `Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]) = Pr(C[i]=c[i] | K[i]=j', T[i]=t[i])`
#'     for any `j, j' in c[i]`.
#' C3: masking probabilities are independent of `theta`
#'
#' As a special case, this model also includes exact component cause of failure
#' data where the candidate set is a singleton.
#' 
#' @author Alex Towell
#' @name Exponential series likelihood model for masked component cause of
#' failure (candidate sets) satisfying conditions C1, C2, and C3 with exact
#' system failure time observations.
#' @keywords exponential, distribution, series, statistics, masked data
exp_series_md_c1_c3 <- function(rates = NULL) {
    structure(list(rates = rates),
              class = c("exp_series_md_c1_c3",
                        "series_md",
                        "likelihood_model"))
}

#' Generates a log-likelihood function for an exponential series system with
#' respect to rate parameter for masked data with candidate sets that satisfy 
#' conditions C1 and C3.
#'
#' @return log-likelihood function that takes the following arguments:
#'  - `df`: masked data
#'  - `theta`: rate / scale parameters of component lifetime distributions
#'  - `sys.var`: system lifetime (optionally right-censored) column name,
#'               defaults to `t`
#'  - `candset`: prefix of Boolean matrix encoding of candidate sets, defaults
#'               to `x`, e.g., `x1,...,xm`.
#' @importFrom md.tools md_decode_matrix
#' @export
loglik.exp_series_md_c1_c3 <- function(model) {
    
    function(df, theta, sys.var = "t", cand.var = "x", ...) {
        if (any(theta <= 0)) return(NA)
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }
        
        f <- 0
        stopifnot(sys.var %in% colnames(df))
        sum.t <- sum(df[[sys.var]])
        C <- md_decode_matrix(df, candset)
        m <- ncol(C)
        if (m == 0) {
            stop("No candidate sets wih prefix '", cand.var, "' found")
        }
        f <- -sum.t * sum(theta)
        for (i in seq_len(n)) {
            f <- f + log(sum(theta[C[i, ]]))
        }
        return(f)
    }
}

#' Generates a score function for an exponential series system (or
#' competing risks) with respect to parameter `theta` for masked component
#' failure (or masked competing risks) with candidate sets (risks) that
#' approximately satisfy conditions C1 and C3.
#'
#' @param df failure-time data with masked competing risks
#' @param P masking probability P{C[i] | K[i], T[i]}
#' @return score function
#' @importFrom md.tools md_decode_matrix
#' @export
score.exp_series_md_c1_c3 <- function(df, P, ...) {

    n <- nrow(df)
    if (n == 0) {
        stop("df is empty")
    }

    function(theta, sys.var = "t", candset = "x") {
        if (length(theta) != m) stop("length(theta) != m")
        if (any(theta <= 0)) return(NA)
        for (j in seq_len(m)) {
            for (i in seq_len(n)) {
                c <- C[i,]
                if (c[j]) {
                    denom <- 0
                    for (r in seq_len(m)) {
                        if (c[r]) {
                            denom <- denom + theta[r] * P(c, r, t[i])
                        }
                    }
                    v[j] <- v[j] + P(c, j, t[i]) / denom
                }
            }
        }
        v
    }
}
