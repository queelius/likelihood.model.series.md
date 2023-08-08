#' Likelihood contribution for exponential series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3,
#' desribed below. 
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
        stop("no candidate sets wih prefix '", cand.var, "' found")
    }
    list(sum.t = sum.t, C = df$C)
}

#' Generates a log-likelihood function for an exponential series system with
#' respect to rate parameter for masked data with candidate sets that satisfy 
#' conditions C1, C2, and C3.
#'
#' @param md masked data
#' @param options list of options
#' - `set.var` prefix of Boolean matrix encoding of candidate sets, defaults
#'   to `x`, e.g., `x1,...,xm`.
#' - `sys.var` system lifetime (optionally right-censored) column name,
#'   defaults to `t`
#' - `delta.var` right-censoring indicator column name. If `NULL`, then
#'   no right-censoring is assumed. If a system lifetime is
#'   right-censored (*not* observed), the right-censoring
#'   indicator is `TRUE`, otherwise it is `FALSE`.
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
            m <- ncol(df$C)
            stopifnot(length(theta) == m)
            f <- -df$sum.t * sum(theta)
            for (i in seq_len(n)) {
                f <- f + df$n[i] * log(sum(theta[df$C[i, ]]))
            }
        } else {
            stopifnot(sys.var %in% colnames(df))
            sum.t <- sum(df[[sys.var]])
            C <- md_decode_matrix(df, cand.var)
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
#' - `set.var` prefix of Boolean matrix encoding of candidate sets, defaults
#'   to `x`, e.g., `x1,...,xm`.
#' - `sys.var` system lifetime (optionally right-censored) column name,
#'   defaults to `t`
#' - `delta.var` right-censoring indicator column name. If `NULL`, then
#'   no right-censoring is assumed. If a system lifetime is
#'   right-censored (*not* observed), the right-censoring
#'   indicator is `TRUE`, otherwise it is `FALSE`.
#' @param ... pass additional arguments to `options`
#' @return score function of type `R^m -> R`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
score.exp_series_system_md_C1_C2_C3 <- function(model) {
    
    function(df, theta, sys.var = "t", cand.var = "x", use_mss = FALSE, ...) {
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }
        stopifnot(sys.var %in% colnames(df))

        if (use_mss) {
            sum.t <- sum(md[ , options$sys.var])
            if (!is.null(options$delta.var) && options$delta.var %in% colnames(md)) {
                # only keep the observations that were not right-censored in `md`
                md <- md %>% filter(!.[[options$delta.var]])
            }
            md$C <- md_decode_matrix(md, options$set.var)
            md <- md %>% group_by(.data$C) %>% count()
            m <- ncol(md$C)
            if (m == 0) {
                stop("no candidate sets wih prefix '", options$set.var, "' found")
            }

            function(theta) {
                if (length(theta) != m) stop("length(theta) != m")
                if (any(theta <= 0)) return(NA)
                v <- rep(-sum.t, m)
                for (j in seq_len(m)) {
                    for (i in seq_len(n)) {
                        if (md$C[i, j]) {
                            v[j] <- v[j] + md$n[i] / sum(theta[md$C[i, ]])
                        }
                    }
                }
                v
            }
        }
    }
}

#' Generates the observed information matrix (FIM) for an exponential series
#' system with respect to parameter `theta` for masked data with candidate
#' sets that approximately satisfy conditions C1, C2, and C3 with exact
#' failure times, or right-censoring (no component failure).
#'
#' @param md masked data with candidate sets that meet the
#'           regular candidate model
#' @return observed information matrix of type `R^m -> R^(m x m)`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
hess_loglik.exp_series_system_md_C1_C2_C3 <- function(model) {

    
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    defaults <- list(
        set.var = "x",
        sys.var = "t",
        delta.var = NULL)

    options <- modifyList(defaults, options)
    options <- modifyList(options, list(...))
    stopifnot(
        options$sys.var %in% colnames(md),
        is.null(options$delta.var) || options$delta.var %in% colnames(md))

    if (!is.null(options$delta.var) && options$delta.var %in% colnames(md))
        md <- md %>% filter(.data$delta == FALSE)
    md$C <- md_decode_matrix(md, options$set.var)
    md <- md %>% group_by(.data$C) %>% count()
    m <- ncol(md$C)
    if (m == 0) {
        stop("no candidate sets wih prefix '", options$set.var, "' found")
    }

    function(df = NULL, theta) {
        if (is.null(df)) {
            df <- md
        } else {
            
        }
        if (length(theta) != m) stop("length(theta) != m")
        if (any(theta <= 0)) return(NA)
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
        I
    }
}
