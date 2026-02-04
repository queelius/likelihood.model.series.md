#' Constructs a likelihood model for `wei_series_homogeneous_md_c1_c2_c3`.
#'
#' Likelihood model for Weibull series systems with homogeneous shape parameter
#' and masked component cause of failure with candidate sets that satisfy
#' conditions C1, C2, and C3.
#'
#' This is a reduced model where all components share a common shape parameter k,
#' while retaining individual scale parameters. The parameter vector is
#' (k, scale_1, ..., scale_m), giving m+1 parameters instead of 2m.
#'
#' A key property of this model is that the series system lifetime is itself
#' Weibull distributed with shape k and scale lambda_s = (sum(scale_j^{-k}))^{-1/k}.
#'
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood.model` package by providing the following methods:
#'
#'  (1) `loglik.wei_series_homogeneous_md_c1_c2_c3`
#'  (2) `score.wei_series_homogeneous_md_c1_c2_c3`
#'  (3) `hess_loglik.wei_series_homogeneous_md_c1_c2_c3`
#'
#' In this likelihood model, masked component data approximately satisfies:
#'
#' C1: `Pr{K[i] in C[i]} = 1`
#' C2: `Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]} = Pr{C[i]=c[i] | K[i]=j', T[i]=t[i]}`
#'     for any `j, j' in c[i]`.
#' C3: masking probabilities are independent of `theta`
#'
#' @param shape common shape parameter for all Weibull component lifetimes (optional)
#' @param scales scale parameters for Weibull component lifetimes (optional)
#' @param lifetime column name for system lifetime, defaults to `"t"`
#' @param indicator column name for right-censoring indicator, defaults to `"delta"`.
#'        TRUE/1 = exact failure time, FALSE/0 = right-censored.
#'        For backwards compatibility, if this column is not present in the data,
#'        censoring is inferred from empty candidate sets (all FALSE).
#' @param candset column prefix for candidate set indicators, defaults to `"x"`
#' @export
#' @return likelihood model object with class
#'         `c("wei_series_homogeneous_md_c1_c2_c3", "series_md", "likelihood_model")`
#' @seealso [wei_series_md_c1_c2_c3()] for the full model with heterogeneous shapes
#' @examples
#' # Create model and fit to data using generic dispatch
#' model <- wei_series_homogeneous_md_c1_c2_c3()
#' # solver <- fit(model)
#' # theta: (shape, scale1, scale2, ...)
#' # mle <- solver(data, par = c(1.2, 1000, 900, 850))
wei_series_homogeneous_md_c1_c2_c3 <- function(shape = NULL, scales = NULL,
                                                lifetime = "t",
                                                indicator = "delta",
                                                candset = "x") {
    structure(
        list(
            shape = shape,
            scales = scales,
            lifetime = lifetime,
            indicator = indicator,
            candset = candset
        ),
        class = c("wei_series_homogeneous_md_c1_c2_c3",
                  "series_md",
                  "likelihood_model")
    )
}


#' Log-likelihood method for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns a log-likelihood function for a Weibull series system with
#' homogeneous shape parameter. The parameter vector is (k, scale_1, ..., scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return log-likelihood function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: parameter vector (shape, scale1, scale2, ...)
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `indicator`: right-censoring indicator column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model loglik
#' @method loglik wei_series_homogeneous_md_c1_c2_c3
#' @export
loglik.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_indicator <- model$indicator %||% "delta"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, indicator = default_indicator,
             candset = default_candset, ...) {
        d <- extract_model_data(df, lifetime, indicator, candset)

        p <- length(par)
        if (p != d$m + 1) stop("length(par) must equal number of components + 1")

        k <- par[1]          # common shape
        scales <- par[-1]    # individual scales

        if (k <= 0 || any(scales <= 0)) return(-Inf)

        # Log-likelihood for Weibull series with homogeneous shape (C1, C2, C3)
        # S(t) = exp(-sum((t/scale_j)^k)),  h_j(t) = (k/scale_j) * (t/scale_j)^(k-1)
        s <- 0
        for (i in seq_len(d$n)) {
            s <- s - sum((d$t[i] / scales)^k)

            if (d$delta[i]) {
                cind <- d$C[i, ]
                if (any(cind)) {
                    hazard_sum <- sum(k / scales[cind] *
                                      (d$t[i] / scales[cind])^(k - 1))
                    if (hazard_sum > 0) {
                        s <- s + log(hazard_sum)
                    }
                }
            }
        }
        return(s)
    }
}


#' Score method for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns a score (gradient) function for a Weibull series system with
#' homogeneous shape parameter. The parameter vector is (k, scale_1, ..., scale_m)
#' for masked data with candidate sets that satisfy conditions C1, C2, and C3.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return score function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: parameter vector (shape, scale1, scale2, ...)
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `indicator`: right-censoring indicator column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model score
#' @method score wei_series_homogeneous_md_c1_c2_c3
#' @export
score.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_indicator <- model$indicator %||% "delta"
    default_candset <- model$candset %||% "x"

    function(df, par, lifetime = default_lifetime, indicator = default_indicator,
             candset = default_candset, ...) {
        d <- extract_model_data(df, lifetime, indicator, candset)

        p <- length(par)
        if (p != d$m + 1) stop("length(par) must equal number of components + 1")

        k <- par[1]          # common shape
        scales <- par[-1]    # individual scales

        if (k <= 0 || any(scales <= 0)) return(rep(NA, p))

        shape_score <- 0
        scale_scores <- rep(0, d$m)

        for (i in seq_len(d$n)) {
            # Let u_j = t/scale_j
            u <- d$t[i] / scales
            u_k <- u^k

            # Survival contribution to score
            # d/dk [-sum(u_j^k)] = -sum(u_j^k * log(u_j))
            shape_score <- shape_score - sum(u_k * log(u))

            # d/d(scale_j) [-sum(u_j^k)] = k * u_j^k / scale_j
            scale_scores <- scale_scores + k * u_k / scales

            # Candidate set contribution (only for exact observations)
            if (d$delta[i]) {
                cind <- d$C[i, ]
                if (any(cind)) {
                    u_km1 <- u[cind]^(k - 1)
                    h_j <- k / scales[cind] * u_km1
                    sum_h <- sum(h_j)

                    if (sum_h > 0) {
                        # d(h_j)/dk = (1/scale_j) * u_j^(k-1) * (1 + k * log(u_j))
                        dh_dk <- (1 / scales[cind]) * u_km1 * (1 + k * log(u[cind]))
                        shape_score <- shape_score + sum(dh_dk) / sum_h

                        # d(h_j)/d(scale_j) = -k * h_j / scale_j
                        dh_dscale <- -k * h_j / scales[cind]
                        scale_scores[cind] <- scale_scores[cind] + dh_dscale / sum_h
                    }
                }
            }
        }

        return(c(shape_score, scale_scores))
    }
}


#' Hessian of log-likelihood method for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns the Hessian (second derivative matrix) of the log-likelihood for a
#' Weibull series system with homogeneous shape. Computed numerically via the
#' Jacobian of the score.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return hessian function that takes the following arguments:
#'  - `df`: masked data frame
#'  - `par`: parameter vector (shape, scale1, scale2, ...)
#'  - `lifetime`: system lifetime column name (default from model)
#'  - `indicator`: right-censoring indicator column name (default from model)
#'  - `candset`: prefix of Boolean matrix encoding candidate sets (default from model)
#' @importFrom numDeriv jacobian
#' @importFrom likelihood.model hess_loglik
#' @method hess_loglik wei_series_homogeneous_md_c1_c2_c3
#' @export
hess_loglik.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
    # Get the score function
    score_fn <- score.wei_series_homogeneous_md_c1_c2_c3(model, ...)

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


#' Assumptions for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns the assumptions made by this likelihood model.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (ignored)
#' @return character vector of assumptions
#' @importFrom likelihood.model assumptions
#' @method assumptions wei_series_homogeneous_md_c1_c2_c3
#' @export
assumptions.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
    c(
        "iid observations",
        "Weibull component lifetimes with COMMON shape parameter",
        "series system configuration",
        "system lifetime is Weibull with shape k and scale (sum(scale_j^{-k}))^{-1/k}",
        "C1: failed component is in candidate set with probability 1",
        "C2: uniform probability for candidate sets given component cause",
        "C3: masking probabilities independent of system parameters"
    )
}


#' System scale parameter for homogeneous Weibull series
#'
#' For a series system with Weibull components sharing shape k but with
#' individual scales, the system lifetime is itself Weibull with shape k
#' and this computed scale.
#'
#' @param k common shape parameter
#' @param scales vector of component scale parameters
#' @return system scale parameter
#' @export
#' @examples
#' # 3-component system with common shape 1.2
#' wei_series_system_scale(k = 1.2, scales = c(1000, 900, 850))
wei_series_system_scale <- function(k, scales) {
    (sum(scales^(-k)))^(-1/k)
}


#' Random data generation for `wei_series_homogeneous_md_c1_c2_c3` model.
#'
#' Returns a function that generates random masked data from the homogeneous
#' Weibull series system DGP at a given parameter value.
#'
#' @param model the likelihood model object
#' @param ... additional arguments (passed to returned function)
#' @return function that takes (theta, n, tau, p, ...) and returns a data frame
#'         with columns: t, delta, x1, x2, ..., xm
#' @importFrom likelihood.model rdata
#' @importFrom stats rweibull runif
#' @method rdata wei_series_homogeneous_md_c1_c2_c3
#' @export
rdata.wei_series_homogeneous_md_c1_c2_c3 <- function(model, ...) {
    # Capture model defaults
    default_lifetime <- model$lifetime %||% "t"
    default_indicator <- model$indicator %||% "delta"
    default_candset <- model$candset %||% "x"

    function(theta, n, tau = Inf, p = 0, ...) {
        if (length(theta) < 2) stop("theta must have at least 2 elements (k, scale1, ...)")

        k <- theta[1]          # common shape
        scales <- theta[-1]    # individual scales
        m <- length(scales)

        if (k <= 0 || any(scales <= 0)) {
            stop("Shape and all scale parameters must be positive")
        }

        # Generate component lifetimes (all with same shape k)
        comp_lifetimes <- matrix(nrow = n, ncol = m)
        for (j in seq_len(m)) {
            comp_lifetimes[, j] <- rweibull(n, shape = k, scale = scales[j])
        }

        generate_masked_series_data(comp_lifetimes, n, m, tau, p,
                                     default_lifetime, default_indicator,
                                     default_candset)
    }
}
