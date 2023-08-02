#' Likelihood model for exponential series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3,
#' desribed below. It also includes right-censoring, but this is not the focus
#' of the model and normally we will disentangle the two and model them
#' as separate and independent contributions to the total likelihood model
#' under consideration.
#' 
#' Various functions are provided according to this model:
#' 
#'  (1) `loglik.exp_series_system_md_C1_C2_C3`
#'  (2) `score.exp_series_system_md_C1_C2_C3`
#'  (3) `hess_loglik.exp_system_series_md_C1_C2_C3`
#'  (4) `fim.exp_series_system_md_C1_C2_C3`
#' 
#' These are useful for doing maximum likelihood estimation, hypothesis
#' testing (e.g., likelihood ratio test), estimation of asymptotic sampling
#' distribution given data from the DGP according to the specified model,
#' etc.
#' 
#' It is designed to work well with the `algebraic.mle` package.
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
#' In this likelihood model, we also assume that the failure time is exactly
#' known.
#' 
#' @author Alex Towell
#' @name Exponential series likelihood model for masked component cause of
#' failure (candidate sets) satisfying conditions C1, C2, and C3 with exact
#' system failure time observations.
#' @keywords exponential, distribution, series, statistics, masked data
exp_series_system_md_C1_C2_C3 <- function(rates = NULL) {
    structure(list(rates = rates),
              class = c("exp_series_system_md_C1_C2_C3",
                        "likelihood_model"))
}

#' Generates a set of minimally sufficient statistic for an exponential series
#' system with respect to rate parameter for masked data with candidate sets
#' that satisfy conditions C1, C2, and C3.
#' 
#' @param df data frame
#' @param sys.var system lifetime column name, defaults to "t"
#' @param cand.var prefix of Boolean matrix encoding of candidate sets, defaults
#' to `x`, e.g., `x1,...,xm`.
#' @return sufficient statistics
exp_series_system_md_C1_C2_C3_mss <- function(df, sys.var = "t", cand.var = "x") {
    n <- nrow(df)
    if (n == 0) {
        stop("df is empty")
    }

    stopifnot(sys.var %in% colnames(df))
    sum.t <- sum(df[[sys.var]])
    df$C <- md_decode_matrix(df, cand.var)
    df <- df %>% group_by(.data$C) %>% count()
    m <- ncol(df$C)
    if (m == 0) {
        stop("no candidate sets found")
    }
    list(sum.t = sum.t, C = df$C, n = df$n)
}

#' Generates a log-likelihood function for an exponential series system with
#' respect to rate parameter for masked data with candidate sets that satisfy 
#' conditions C1, C2, and C3.
#'
#' @param md masked data
#' @param options list of options
#' - `candset` prefix of Boolean matrix encoding of candidate sets, defaults
#'   to `x`, e.g., `x1,...,xm`.
#' - `sys.var` system lifetime (optionally right-censored) column name,
#'   defaults to `t`
#' @return log-likelihood function
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
loglik.exp_series_system_md_C1_C2_C3 <- function(model) {
    
    function(df, theta, sys.var = "t", cand.var = "x", use_mss = FALSE, ...) {
        if (any(theta <= 0)) return(NA)
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }
        

        f <- 0
        if (use_mss) {
            stopifnot("sum.t" %in% names(df), "C" %in% names(df), "n" %in% names(df))
            m <- ncol(df$C)
            stopifnot(length(theta) == m)
            
            f <- -df$sum.t * sum(theta)
            for (i in seq_len(n)) {
                f <- f + df$n[i] * log(sum(theta[df$C[i, ]]))
            }
        } else {
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
        }
        return(f)
    }
}


#' Generates a score function for an exponential series system (or
#' competing risks) with respect to parameter `theta` for masked component
#' failure (or masked competing risks) with candidate sets (risks) that
#' approximately satisfy conditions C1, C2, and C3.
#'
#' @param md right-censored failure-time data with masked competing risks
#' @param options list of options
#' - `candset` prefix of Boolean matrix encoding of candidate sets, defaults
#'   to `x`, e.g., `x1,...,xm`.
#' - `sys.var` system lifetime (optionally right-censored) column name,
#'   defaults to `t`
#' @param ... pass additional arguments to `options`
#' @return score function of type `R^m -> R`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
score.exp_series_system_md_C1_C2_C3 <- function(model) {
    
    function(df, theta, sys.var = "t", candset = "x", use_mss = FALSE, ...) {
        if (any(theta <= 0)) return(NA)
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }
        stopifnot(sys.var %in% colnames(df))

        if (use_mss) {
            stopifnot("sum.t" %in% names(df), "C" %in% names(df), "n" %in% names(df))
            m <- ncol(df$C)
            if (length(theta) != m) stop("length(theta) != m")        

            v <- rep(-df$sum.t, m)
            for (j in seq_len(m)) {
                for (i in seq_len(n)) {
                    if (df$C[i, j]) {
                        v[j] <- v[j] + df$n[i] / sum(theta[df$C[i, ]])
                    }
                }
            }
        } else {
            stopifnot(sys.var %in% colnames(df), candset %in% colnames(df))
            C <- md_decode_matrix(df, candset)
            m <- ncol(C)
            stopifnot(length(theta) == m)

            v <- rep(-sum(df[[sys.var]]), m)
            for (j in seq_len(m)) {
                for (i in seq_len(n)) {
                    if (md$C[i, j]) {
                        v[j] <- v[j] + md$n[i] / sum(theta[md$C[i, ]])
                    }
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
#' @param md masked data with candidate sets that meet the
#'           regular candidate model
#' @return hessian of the loglikelihood of `loglik.exp_series_system_md_c1_c2_c3`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
hess_loglik.exp_series_system_md_C1_C2_C3 <- function(model) {

    function(df, theta, sys.var = "t", candset = "x" ...) {

        if (any(theta <= 0)) return(NA)

        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }

        stopifnot(sys.var %in% colnames(df))

        if (use_mss) {
        
        } else {
            df$C <- md_decode_matrix(df, candset)
            m <- ncol(df$C)
            if (m == 0) {
                stop("no candidate sets found")
            }
            if (length(theta) != m) stop("length(theta) != m")

            df <- df %>% group_by(.data$C) %>% count()
            I <- matrix(rep(0, m * m), nrow = m)
            for (j in 1:m) {
                for (k in 1:m) {
                    for (i in 1:n) {
                        if (md$C[i, j] && md$C[i, k]) {
                            I[j, k] <- I[j, k] + md$n[i] / sum(theta[md$C[i, ]])^2
                        }
                    }
                }
            }
        }
    }
}