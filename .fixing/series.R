#' This file contains functions related to a general series distribution.
#' Functions include simulation, pdf, cdf, quantile, and other related
#' functions. It models the concept of an Algebraic Distribution, in that
#' implements some generic methods in `algebraic.dist`.
#' 
#' It is based on the paper:
#'
#' Reliability Estimation in Series Systems: Maximum Likelihood Techniques for Right-Censored and Masked Failure Data
#' GitHub link: \url{https://github.com/queelius/reliability-estimation-in-series-systems}.
#' 
#' For a likelihood model (for parameter estimation) from masked data
#' (competing risks) from candidate sets satisfying conditions C1, C2,
#' and C3, see the file `series_md_c1_c2_c3.R`.
#'
#' @author Alex Towell
#' @name General series system
#' @keywords distribution, series, statistics
NULL

#' General series system object constructor.
#'
#' @param hazs list of hazard functions.
#' @param nparam integer vector, we have m components, each of which
#'               may have a different number of parameters. the j-th
#'               component has `nparams[j]` parameters.
#' @param survs list of component survival (reliability) functions.
#'              Defaults to NULL, in which case we generate the survival
#'              functions from `hazs` (using cumulative hazard functions).
#' @return `series` system object
#' @export
series <- function(
    hazs,
    nparams,
    survs = NULL) {

    stopifnot(length(hazs) == length(nparams))

    # use cumulative hazard functions to generate survival functions
    if (is.null(survs)) {
        survs <- lapply(hazs, cum_haz)
    }

    structure(
        list(hazs = hazs,
             survs = survs,
             nparams = nparams),
        class = c("series",
                  "latent_multivariate_dist", 
                  "univariate_dist",
                  "dist"))
}

#' Determine if an object is a `series` object.
#' @param x object to test
#' @return logical, TRUE if `x` is a `series` object
#' @export
is_series <- function(x) {
    inherits(x, "series")
}

#' Hazard function for a general series system.
#' 
#' @param object series object
#' @return hazard function for series system
#' @export
hazard.series <- function(object) {
    h <- hazard_series_helper(object$hazs,
                              object$nparams)
    function(t, theta, ...) {
        h(t, theta, ...)
    }
}

#' Survival function for a general series system.
#' 
#' @param object `series` object
#' @return survival function
#' @importFrom algebraic.dist surv
#' @export
surv.series <- function(object) {
    S <- survival_series_helper(object$survs,
                                object$nparams)

    function(t, theta = NULL, log = FALSE, ...) {
        if (is.null(theta)) {
            theta <- param(object)
        }
        S(t, theta, log, ...)
    }
}


#' Cumulative distribution function (cdf) for a general series system.
#' @param object series object
#' @importFrom algebraic.dist surv cdf
#' @export
cdf.series <- function(object, ...) {
    R <- surv(object, ...)
    function(t, theta, ...) {
        1 - R(t, theta, ...)
    }
}

#' Component hazard function for a general series system.
#' @param x series object
#' @return list of component hazard functions for series system
#' @export
hazards.series <- function(x) {
    x$hazs
}

#' Number of parameters for a general series system.
#' @param x series object
#' @return number of parameters
#' @importFrom algebraic.dist nparams
#' @export
nparams.series <- function(x) {
    sum(x$nparams)
}

#' Quantile function (inverse of the cdf) for a general series system.
#' By definition, the quantile `p * 100%` is the value `t` that
#' satisfies `F(t) - p = 0`. We solve for `t` using Newton's method.
#'
#' @param p vector of probabilities.
#' @param theta parameter vector
#' @param object series object
#' @param ... additional arguments to pass into `options`
#' @importFrom algebraic.dist params inv_cdf
#' @importFrom stats optim
#' @export
inv_cdf.series <- function(x, ...) {

    h <- hazard(x)
    R <- surv(x)
    p <- sum(nparams(x))

    function(p, theta, t0 = 1) {
        stopifnot(length(theta) == p, all(p > 0))
        optim(
            par = t0,
            fn = function(t) R(t, theta) - p,
            gr = function(t) h(t, theta) * R(t, theta),
            method = "Brent",
            lower = 0,
            upper = Inf)$par
    }
}

#' Sample from a general series system.
#'
#' @param n sample size
#' @param theta parameter vector, if NULL uses the values in `object`
#' @param object a `series` object
#' @param options list of options
#' @param ... additional arguments to pass into `options`
#' @export
sampler.series <- function(x, ...) {

    # generate sampling function from survival functions
    # note that T[i] ~ min(T1[i],...,Tm[i]) where Tj[i] ~ Sj
    # so, we sample each component and take the minimum
    # to sample a component, we solve for x in
    # F(x) = p => S(x) - (1-p) = 0
    # where p ~ UNIF(0,1) using Newton's method.
    for (j in 1:m) {
        # theta_comp[j] <- theta[cur_index:(cur_index + nparams[j] - 1)]
        samps[[j]] <- function(n) {
            p <- runif(n)
            for (i in 1:n) {
                stats::optim(
                    par = options$t0,
                    fn = function(t) x$survs[[j]](t, theta) - p[i],
                    gr = function(t) x$hazs[[j]](t, theta) / x$survs[[j]](t, theta),
                    method = "Brent",
                    lower = 0,
                    upper = Inf)
            }
        }
    }

    function(n, theta = NULL, options, ...) {

        defaults <- list(
            include_latent = FALSE,
            ...)
        options <- modifyList(defaults, options)

        # Create an empty matrix to store the lifetimes
        lifetimes <- matrix(nrow = n, ncol = length(comps))
        
        # Initialize a counter for the current index in theta
        cur_index <- 1
        # For each component in the series

        m <- ncomponents(x)

        for (i in 1:m) {
            comp <- comps[[i]]
            params <- theta[cur_index:(cur_index + nparams[i] - 1)]
            
            # Update the current index for the next iteration
            cur_index <- cur_index + num_params
            
            # Generate random lifetimes for this component
            lifetimes[,i] <- do.call(sampler, c(list(n), as.list(params)))
        }
        
        data <- apply(lifetimes, 1, min)
        if (options$include_latent) {
            data <- data %>% md_encode_matrix("t") %>%
                            md_mark_latent(paste0("t",1:m))
        }
        data
    }
}

#' Probability density function (pdf) for a general series system.
#' @param object series object
#' @return pdf function
#' @importFrom algebraic.dist surv hazard
#' @importFrom stats density
#' @export
density.series <- function(object) {
    S <- surv(object)
    h <- hazard(object)

    function(t, theta, log.p = FALSE, ...) {
        A <- S(t, theta, log.p = log.p)
        B <- h(t, theta, log.p = log.p)
        ifelse(log.p, A + B, A * B)
    }
}




#' Cumulative distribution function (cdf) for `series` systems.
#'
#' @param x series object
#' 
#' @return cdf function, a function of `t` of times, parameter vector
#' `theta`, and logical `log` indicating whether to return the log of
#' the cdf
#' @export
cdf.series <- function(x) {

    R <- surv(x)
    p <- sum(nparams(x))

    function (t, theta, lower.limit = TRUE, log = FALSE) {

        stopifnot(length(theta) == p)

        if (lower.limit) {
            val <- R(t, theta, log = FALSE)
            if (log)
                return(log(1 - val))
            else
                return(1 - val)
        } else {
            return(R(t, theta, log = log))
        }
    }
}

#' Survival function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param object series object
#' @param log logical, log of survival function
#' @export
surv.series <- function(x) {

    R <- surv(x)
    p <- sum(nparams(x))

    # survival_general_series_helper <- function(
    #     survivals,
    #     param_counts) {

    #     m <- length(survivals)
    #     function(t, theta, log.p = FALSE) {
    #         i0 <- 1
    #         i1 <- param_counts[1]
    #         p <- ifelse(log.p, 0, 1)
    #         for (j in 1:m) {
    #             P <- survivals[[j]](t, theta[i0:i1], log.p = log.p)
    #             p <- ifelse(log.p, p + P, p * P)
    #             i0 <- i1 + 1
    #             i1 <- i1 + param_counts[j]
    #         }
    #         p
    #     }
    # }

    function(t, theta, log = FALSE) {

        stopifnot(length(theta) == p)

        R(t, theta, log = log)
    }
}

#' Hazard function for `series` system
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @export
hazard.series <- function(object) {

    h <- hazard(object)
    p <- sum(nparams(object))

    # hazard_general_series_helper <- function(
    #     hazards,
    #     param_counts) {

    #     m <- length(hazards)
    #     function(t, theta, log.p = FALSE) {
    #         i0 <- 1
    #         i1 <- param_counts[1]
    #         v <- 0
    #         for (j in 1:m) {
    #             v <- v + hazards[[j]](t,
    #                 theta[i0:i1],
    #                 log.p = FALSE)
    #             i0 <- i1 + 1
    #             i1 <- i1 + param_counts[j]
    #         }
    #         ifelse(log.p, log(v), v)
    #     }
    # }


    function(t, theta, log = FALSE) {
        stopifnot(length(theta) == p)
        h(t, theta, log = log)
    }
}






















