#' Constructs a likelihood model for `wei_series_md_c1_c2_c3`.
#'
#' Likelihood model for Weibull series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood.model` package by providing the following methods:
#'
#'  (1) `loglik.wei_series_md_c1_c2_c3`
#'  (2) `score.wei_series_md_c1_c2_c3`
#'  (3) `hess_loglik.wei_series_md_c1_c2_c3`
#'
#' The Weibull series system has 2m parameters: (shape_1, scale_1, ..., shape_m, scale_m).
#'
#' In this likelihood model, masked component data approximately satisfies:
#'
#' C1: `Pr{K[i] in C[i]} = 1`
#' C2: `Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]} = Pr{C[i]=c[i] | K[i]=j', T[i]=t[i]}`
#'     for any `j, j' in c[i]`.
#' C3: masking probabilities are independent of `theta`
#'
#' @param shapes shape parameters for Weibull component lifetimes (optional)
#' @param scales scale parameters for Weibull component lifetimes (optional)
#' @param lifetime column name for system lifetime, defaults to `"t"`
#' @param indicator column name for right-censoring indicator, defaults to `"delta"`.
#'        TRUE/1 = exact failure time, FALSE/0 = right-censored.
#'        For backwards compatibility, if this column is not present in the data,
#'        censoring is inferred from empty candidate sets (all FALSE).
#' @param candset column prefix for candidate set indicators, defaults to `"x"`
#' @export
#' @return likelihood model object with class
#'         `c("wei_series_md_c1_c2_c3", "series_md", "likelihood_model")`
#' @examples
#' # Create model and fit to data using generic dispatch
#' model <- wei_series_md_c1_c2_c3()
#' # solver <- fit(model)
#' # theta: (shape1, scale1, shape2, scale2, ...)
#' # mle <- solver(data, par = c(1, 1000, 1, 1000, 1, 1000))
wei_series_md_c1_c2_c3 <- function(shapes = NULL, scales = NULL, lifetime = "t",
                                   indicator = "delta", candset = "x") {
    structure(
        list(
            shapes = shapes,
            scales = scales,
            lifetime = lifetime,
            indicator = indicator,
            candset = candset
        ),
        class = c("wei_series_md_c1_c2_c3",
                  "series_md",
                  "likelihood_model")
    )
}


#' Log-likelihood method for `wei_series_md_c1_c2_c3` model.
#'
#' Returns a log-likelihood function for a Weibull series system with
#' respect to parameter vector (shape_1, scale_1, ..., shape_m, scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `indicator`: right-censoring indicator column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model loglik
#' @method loglik wei_series_md_c1_c2_c3
#' @export
loglik.wei_series_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_indicator <- model$indicator %||% "delta"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, indicator = default_indicator,
             candset = default_candset, ...) {
        n <- nrow(df)
        if (n == 0) stop("df is empty")
        if (!lifetime %in% colnames(df)) stop("lifetime variable not in colnames(df)")

        C <- md_decode_matrix(df, candset)
        if (is.null(C)) stop("no candidate set found for candset")
        m <- ncol(C)

        k <- length(par)
        if (k != 2 * m) stop("length(par) must equal 2 * number of components")

        shapes <- par[seq(1, k, 2)]
        scales <- par[seq(2, k, 2)]

        if (any(shapes <= 0) || any(scales <= 0)) return(-Inf)

        # Get censoring indicator (backwards compatibility)
        if (indicator %in% colnames(df)) {
            delta <- as.logical(df[[indicator]])
        } else {
            delta <- rowSums(C) > 0
        }

        t <- df[[lifetime]]

        # Log-likelihood for Weibull series with masked data (C1, C2, C3)
        s <- 0
        for (i in seq_len(n)) {
            # Survival contribution: -sum((t/scale)^shape)
            s <- s - sum((t[i] / scales)^shapes)

            # Candidate set contribution (only for exact observations)
            if (delta[i]) {
                cind <- C[i, ]
                if (any(cind)) {
                    hazard_sum <- sum(shapes[cind] / scales[cind] *
                                      (t[i] / scales[cind])^(shapes[cind] - 1))
                    if (hazard_sum > 0) {
                        s <- s + log(hazard_sum)
                    }
                }
            }
        }
        return(s)
    }
}


#' Score method for `wei_series_md_c1_c2_c3` model.
#'
#' Returns a score (gradient) function for a Weibull series system with
#' respect to parameter vector (shape_1, scale_1, ..., shape_m, scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return score function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `indicator`: right-censoring indicator column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model score
#' @method score wei_series_md_c1_c2_c3
#' @export
score.wei_series_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_indicator <- model$indicator %||% "delta"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, indicator = default_indicator,
             candset = default_candset, ...) {
        n <- nrow(df)
        if (n == 0) stop("df is empty")
        if (!lifetime %in% colnames(df)) stop("lifetime variable not in colnames(df)")

        C <- md_decode_matrix(df, candset)
        if (is.null(C)) stop("no candidate set found for candset")
        m <- ncol(C)

        k <- length(par)
        if (k != 2 * m) stop("length(par) must equal 2 * number of components")

        shapes <- par[seq(1, k, 2)]
        scales <- par[seq(2, k, 2)]

        if (any(shapes <= 0) || any(scales <= 0)) return(rep(NA, k))

        # Get censoring indicator (backwards compatibility)
        if (indicator %in% colnames(df)) {
            delta <- as.logical(df[[indicator]])
        } else {
            delta <- rowSums(C) > 0
        }

        t <- df[[lifetime]]
        shape_scores <- rep(0, m)
        scale_scores <- rep(0, m)

        for (i in seq_len(n)) {
            # Survival contribution to score
            rt_term_shapes <- -(t[i] / scales)^shapes * log(t[i] / scales)
            rt_term_scales <- (shapes / scales) * (t[i] / scales)^shapes

            # Candidate set contribution (only for exact observations)
            mask_term_shapes <- rep(0, m)
            mask_term_scales <- rep(0, m)

            if (delta[i]) {
                cind <- C[i, ]
                if (any(cind)) {
                    denom <- sum(shapes[cind] / scales[cind] *
                                 (t[i] / scales[cind])^(shapes[cind] - 1))

                    if (denom > 0) {
                        numer_shapes <- 1/t[i] * (t[i] / scales[cind])^shapes[cind] *
                            (1 + shapes[cind] * log(t[i] / scales[cind]))
                        mask_term_shapes[cind] <- numer_shapes / denom

                        numer_scales <- (shapes[cind] / scales[cind])^2 *
                            (t[i] / scales[cind])^(shapes[cind] - 1)
                        mask_term_scales[cind] <- numer_scales / denom
                    }
                }
            }

            shape_scores <- shape_scores + rt_term_shapes + mask_term_shapes
            scale_scores <- scale_scores + rt_term_scales - mask_term_scales
        }

        # Interleave shape and scale scores
        scr <- rep(0, k)
        scr[seq(1, k, 2)] <- shape_scores
        scr[seq(2, k, 2)] <- scale_scores
        return(scr)
    }
}


#' Hessian of log-likelihood method for `wei_series_md_c1_c2_c3` model.
#'
#' Returns the Hessian (second derivative matrix) of the log-likelihood for a
#' Weibull series system. Computed numerically via the Jacobian of the score.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return hessian function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: parameter vector (shape1, scale1, shape2, scale2, ...)
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `indicator`: right-censoring indicator column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom numDeriv jacobian
#' @importFrom likelihood.model hess_loglik
#' @method hess_loglik wei_series_md_c1_c2_c3
#' @export
hess_loglik.wei_series_md_c1_c2_c3 <- function(model, ...) {
    # Get the score function
    score_fn <- score.wei_series_md_c1_c2_c3(model, ...)

    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_indicator <- model$indicator %||% "delta"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, indicator = default_indicator,
             candset = default_candset, ...) {
        numDeriv::jacobian(
            func = function(theta) score_fn(df, theta, lifetime, indicator, candset),
            x = par,
            method.args = list(r = 6)
        )
    }
}


#' Assumptions for `wei_series_md_c1_c2_c3` model.
#'
#' Returns the assumptions made by this likelihood model.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (ignored)
#' @return character vector of assumptions
#' @importFrom likelihood.model assumptions
#' @method assumptions wei_series_md_c1_c2_c3
#' @export
assumptions.wei_series_md_c1_c2_c3 <- function(model, ...) {
    c(
        "iid observations",
        "Weibull component lifetimes",
        "series system configuration",
        "C1: failed component is in candidate set with probability 1",
        "C2: uniform probability for candidate sets given component cause",
        "C3: masking probabilities independent of system parameters"
    )
}
