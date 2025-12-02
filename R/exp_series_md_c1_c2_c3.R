#' Constructs a likelihood model for `exp_series_md_c1_c2_c3`.
#'
#' Likelihood model for exponential series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3,
#' described below.
#'
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood_model` package by providing the following methods:
#'
#'  (1) `loglik.exp_series_md_c1_c2_c3`
#'  (2) `score.exp_series_md_c1_c2_c3`
#'  (3) `hess_loglik.exp_series_md_c1_c2_c3`
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
#' @param rates rate parameters for exponential component lifetimes (optional,
#'              used as initial values for MLE if provided)
#' @param lifetime column name for system lifetime, defaults to `"t"`
#' @param candset column prefix for candidate set indicators, defaults to `"x"`
#' @export
#' @return likelihood model object with class
#'         `c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")`
#' @examples
#' # Create model and fit to data using generic dispatch
#' model <- exp_series_md_c1_c2_c3()
#' # solver <- fit(model)
#' # mle <- solver(data, par = c(1, 1, 1))
exp_series_md_c1_c2_c3 <- function(rates = NULL, lifetime = "t", candset = "x") {
    structure(
        list(
            rates = rates,
            lifetime = lifetime,
            candset = candset
        ),
        class = c("exp_series_md_c1_c2_c3",
                  "series_md",
                  "likelihood_model")
    )
}


#' Log-likelihood method for `exp_series_md_c1_c2_c3` model.
#'
#' Returns a log-likelihood function for an exponential series system with
#' respect to rate parameters for masked data with candidate sets that satisfy
#' conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: rate / scale parameters of component lifetime distributions
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model loglik
#' @method loglik exp_series_md_c1_c2_c3
#' @export
loglik.exp_series_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, candset = default_candset, ...) {
        if (any(par <= 0)) return(-Inf)
        stopifnot(lifetime %in% colnames(df))
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }

        C <- md_decode_matrix(df, candset)
        m <- ncol(C)
        if (m == 0) {
            stop("No candidate sets with prefix '", candset, "' found")
        }

        # Log-likelihood: -sum(t) * sum(theta) + sum_{i: C_i non-empty} log(sum(theta[C_i]))
        # The candidate set term is only included for uncensored observations
        # (right-censored observations have empty candidate sets)
        f <- -sum(df[[lifetime]]) * sum(par)
        for (i in seq_len(n)) {
            theta_c <- sum(par[C[i, ]])
            if (theta_c > 0) {
                f <- f + log(theta_c)
            }
        }
        return(f)
    }
}

#' Score method for `exp_series_md_c1_c2_c3` model.
#'
#' Returns a score (gradient) function for an exponential series system with
#' respect to parameter `theta` for masked component failure with candidate
#' sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return score function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: rate / scale parameters of component lifetime distributions
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model score
#' @method score exp_series_md_c1_c2_c3
#' @export
score.exp_series_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, candset = default_candset, ...) {
        if (any(par <= 0)) return(rep(NA, length(par)))
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }
        stopifnot(lifetime %in% colnames(df))
        C <- md_decode_matrix(df, candset)
        if (is.null(C)) {
            stop("No candidate sets found")
        }
        m <- ncol(C)
        stopifnot(length(par) == m)

        # Score: d/d(theta_j) = -sum(t) + sum(1/sum(theta[C_i]) * I(j in C_i))
        v <- rep(-sum(df[[lifetime]]), m)
        for (j in seq_len(m)) {
            for (i in seq_len(n)) {
                if (C[i, j]) {
                    v[j] <- v[j] + 1 / sum(par[C[i, ]])
                }
            }
        }
        return(v)
    }
}

#' Hessian of log-likelihood method for `exp_series_md_c1_c2_c3` model.
#'
#' Returns the observed information matrix (negative Hessian) for an exponential
#' series system with respect to parameter `theta` for masked data with candidate
#' sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return hessian function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: rate / scale parameters of component lifetime distributions
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model hess_loglik
#' @method hess_loglik exp_series_md_c1_c2_c3
#' @export
hess_loglik.exp_series_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, candset = default_candset, ...) {
        if (any(par <= 0)) return(matrix(NA, length(par), length(par)))
        n <- nrow(df)
        if (n == 0) {
            stop("df is empty")
        }
        C <- md_decode_matrix(df, candset)
        m <- ncol(C)
        if (m == 0) {
            stop("No candidate sets found")
        }
        if (length(par) != m) stop("length(par) != m")

        # Hessian: d^2/d(theta_j)d(theta_k) = -sum(1/sum(theta[C_i])^2 * I(j,k in C_i))
        H <- matrix(0, nrow = m, ncol = m)
        for (j in seq_len(m)) {
            for (k in seq_len(m)) {
                for (i in seq_len(n)) {
                    if (C[i, j] && C[i, k]) {
                        H[j, k] <- H[j, k] - 1 / sum(par[C[i, ]])^2
                    }
                }
            }
        }
        return(H)
    }
}

#' Assumptions for `exp_series_md_c1_c2_c3` model.
#'
#' Returns the assumptions made by this likelihood model.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (ignored)
#' @return character vector of assumptions
#' @importFrom likelihood.model assumptions
#' @method assumptions exp_series_md_c1_c2_c3
#' @export
assumptions.exp_series_md_c1_c2_c3 <- function(model, ...) {
    c(
        "iid observations",
        "exponential component lifetimes",
        "series system configuration",
        "C1: failed component is in candidate set with probability 1",
        "C2: uniform probability for candidate sets given component cause",
        "C3: masking probabilities independent of system parameters"
    )
}

# Null-coalescing operator (if not already defined)
`%||%` <- function(x, y) if (is.null(x)) y else x
