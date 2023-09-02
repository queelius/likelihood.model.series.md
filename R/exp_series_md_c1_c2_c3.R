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
exp_series_md_c1_c2_c3 <- function(rates = NULL) {
    structure(list(rates = rates),
              class = c("exp_series_md_c1_c2_c3",
                        "series_md_c1_c2_c3",
                        "likelihood_model"))
}


#' Generates a log-likelihood function for an exponential series system with
#' respect to rate parameter for masked data with candidate sets that satisfy 
#' conditions C1, C2, and C3.
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
loglik.exp_series_md_c1_c2_c3 <- function(model) {
    
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
            stop("no candidate sets wih prefix '", cand.var, "' found")
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
#' failure (or masked competing risks) with candidate sets (competing risks)
#' that satisfy conditions C1, C2, and C3.
#'
#' @param model the exponential series system masked data model (competing risks)
#'              for candidate set model that satisfies C1, C2, and C3
#' @return score function that takes the following arguments:
#'  - `df`: masked data
#'  - `theta`: rate / scale parameters of component lifetime distributions
#'  - `sys.var`: system lifetime (optionally right-censored) column name,
#'               defaults to `t`
#'  - `candset`: prefix of Boolean matrix encoding of candidate sets, defaults
#'               to `x`, e.g., `x1,...,xm`.
#' @importFrom md.tools md_decode_matrix
#' @export
score.exp_series_md_c1_c2_c3 <- function(model) {
    
    function(df, theta, sys.var = "t", candset = "x") {
        if (any(theta <= 0)) return(NA)
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }
        stopifnot(sys.var %in% colnames(df))
        C <- md_decode_matrix(df, candset)
        if (is.null(C)) {
            stop("no candidate sets found")
        }
        m <- ncol(C)
        stopifnot(length(theta) == m)
        v <- rep(-sum(df[[sys.var]]), m)
        for (j in seq_len(m)) {
            for (i in seq_len(n)) {
                # condition C1: Pr{K[i] in C[i]) = 1
                if (df$C[i, j]) {
                    v[j] <- v[j] + 1 / sum(theta[df$C[i, ]])
                }
            }
        }
        return(v)
    }
}

#' Generates the observed information matrix (FIM) for an exponential series
#' system with respect to parameter `theta` for masked data with candidate
#' sets that approximately satisfy conditions C1, C2, and C3 with exact
#' failure times, or right-censoring (no component failure).
#'
#' @param model the exponential series system masked data model (competing risks)
#'              for candidate set model that satisfies C1, C2, and C3
#' @return hessian of the log-likelihood function `loglik.exp_series_md_c1_c2_c3`
#'         that takes the following arguments:
#'  - `df`: masked data
#'  - `theta`: rate / scale parameters of component lifetime distributions
#'  - `sys.var`: system lifetime (optionally right-censored) column name,
#'               defaults to `t`
#'  - `candset`: prefix of Boolean matrix encoding of candidate sets, defaults
#'               to `x`, e.g., `x1,...,xm`.
#' @importFrom md.tools md_decode_matrix
#' @export
hess_loglik.exp_series_md_C1_C2_C3 <- function(model) {

    function(df, theta, sys.var = "t", candset = "x") {

        if (any(theta <= 0)) return(NA)

        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }

        stopifnot(sys.var %in% colnames(df))

        df$C <- md_decode_matrix(df, candset)
        m <- ncol(df$C)
        if (m == 0) {
            stop("no candidate sets found")
        }
        if (length(theta) != m) stop("length(theta) != m")

        I <- matrix(rep(0, m * m), nrow = m)
        for (j in 1:m) {
            for (k in 1:m) {
                for (i in 1:n) {
                    if (df$C[i, j] && df$C[i, k]) {
                        I[j, k] <- I[j, k] + 1 / sum(theta[df$C[i, ]])^2
                    }
                }
            }
        }
        return(I)
    }
}