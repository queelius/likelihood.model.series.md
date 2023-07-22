#' General series system
#'
#' This file contains functions related to a general series distribution.
#' Functions include simulation, pdf, cdf, quantile, and other related
#' functions.
#' 
#' It is based strictly on the paper "Series system estimation from masked
#' data" by Alex Towell. The paper is available at \url{?}.
#' 
#' For parameter estimation from masked data, see the file
#' `md_series_system.R`.
#'
#' @author Alex Towell
#' @name General series
#' @keywords distribution, series, statistics
NULL


#' General series system object constructor.
#'
#' @param hazs list of hazard functions.
#' @param survs list of reliability functions.
#' @param nparam integer vector, we have m components, each of which
#'                     may have a different number of parameters. the j-th
#'                     component has `nparams[j]` parameters.
#' @return series_system object
#' @export
series_system <- function(
    hazs,
    nparams,
    survs = NULL,
    samps = NULL) {
        stopifnot(length(hazs) == length(nparams))

        if (is.null(survs)) {
            survs <- lapply(hazs, function(h) {
                function(t, theta, ...) {
                    exp(-h(t, theta, ...))
                }
            })
        }

        structure(
            list(hazs = hazs,
                 survs = survs,
                 samps = samps,
                 nparams = nparams),
            class = c("series_system",
                      "latent_multivariate_dist", 
                      "univariate_dist", "dist"))
}

#' Hazard function for a general series system.
#' 
#' @param object series_system object
#' @return hazard function for series system
#' @export
hazard.series_system <- function(object) {
    h <- hazard_series_helper(object$hazs,
                              object$nparams)
    function(t, theta, ...) {
        h(t, theta, ...)
    }
}

#' Survival function for a general series system.
#' 
#' @param object series_system object
#' @return survival function
#' @export
surv.series_system <- function(object) {
    S <- survival_series_helper(object$survs,
                                object$nparams)

    function(t, theta = NULL, log = FALSE, ...) {
        if (is.null(theta)) {
            theta <- param(object)
        }
        S(t, theta, log, ...)
    }
}

#' Probability density function (pdf) for a general series system.
#' @param object series_system object
#' @return pdf function
#' @export
pdf.series_system <- function(object) {
    S <- surv(object)
    h <- hazard(object)

    function(t, theta, log.p = FALSE, ...) {
        A <- S(t, theta, log.p = log.p)
        B <- h(t, theta, log.p = log.p)
        ifelse(log.p, A + B, A * B)
    }
}

#' Cumulative distribution function (cdf) for a general series system.
#' @param object series_system object
#' @export
cdf.series_system <- function(object) {
    R <- surv(object)
    function(t, theta, ...) {
        1 - R(t, theta, ...)
    }
}

#' Quantile function (inverse of the cdf) for a general series system.
#' 
#' By definition, the quantile `p * 100%` is the value `t` that
#' satisfies `F(t) - p = 0`. We solve for `t` using Newton's method.
#' 
#' @param object series_system object
#' @return quantile function, a function of vector `p` of probabilities
#' @export
inv_cdf.series_system <- function(object) {
    Q <- qseries_system(object)
    function(p, theta, ...) {
        Q(p, theta, ...)
    }
}

#' Component hazard function for a general series system.
#' @param object series_system object
#' @return list of component hazard functions for series system
#' @export
hazards.series_system <- function(object) {
    object$hazards
}

#' Number of parameters for a general series system.
#' @param object series_system object
#' @return number of parameters
#' @export
nparams.series_system <- function(object) {
    sum(object$nparams)
}

#' Quantile function (inverse of the cdf) for a general series system.
#' By definition, the quantile `p * 100%` is the value `t` that
#' satisfies `F(t) - p = 0`. We solve for `t` using Newton's method.
#'
#' @param p vector of probabilities.
#' @param theta parameter vector
#' @param object series_system object
#' @param options list of options
#' @param ... additional arguments to pass into `options`
#' @export
inv_cdf.series_system <- function(x, ...) {

    h <- hazard(x)
    R <- surv(x)
    p <- sum(nparams(x))

    function(p, theta = NULL, t0 = 1) {
        theta <- params(x, theta)
        stopifnot(length(theta) == p, all(p > 0))

        optim(
            par = options$t0,
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
#' @param object a `series_system` object
#' @param options list of options
#' @param ... additional arguments to pass into `options`
#' @export
sampler.series_system <- function(x, ...) {

    function(n, theta = NULL, options, ...) {

        defaults <- list(
            include_latent = FALSE,
            ...)
        options <- modifyList(defaults, options)

        stopifnot(is_series_system(object))

        #### previous: qgen_series(runif(n),theta,nparams,h,R)

        comps <- components(object)

        # Create an empty matrix to store the lifetimes
        lifetimes <- matrix(nrow = n, ncol = length(comps))
        
        # Initialize a counter for the current index in theta
        cur_index <- 1
        # For each component in the series
        for (i in seq_along(comps)) {
            comp <- comps[[i]]
            
            # Get the distribution function
            sampler <- match.fun(paste0("r", comp$dist_name))
            
            # Get the parameters for this component
            num_params <- comp$num_params
            params <- theta[cur_index:(cur_index + num_params - 1)]
            
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

#' @export
is_valid_series_system <- function(x) {
    !is.null(x) && is_series_system(x) && !is.null(x$hazards) &&
        length(x$hazards)==length(x$survs)
}

#' pdf for general series system
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param object series_system object
#' @export
pdf.series_system <- function(x) {

    h <- hazard(x)
    R <- surv(x)
    p <- sum(nparams(x))

    function(t, theta = NULL, log = FALSE) {

        theta <- params(x, theta)
        stopifnot(length(theta) == p)

        if (log)
            return(h(t,theta,log=log) - R(t,theta,log=log))
        else
            return(h(t,theta,log=log) / R(t,theta,log=log))

    }
}

#' cdf for general series system
#'
#' @param x series_system object
#' 
#' @return cdf function, a function of `t` of times, parameter vector
#' `theta`, and logical `log` indicating whether to return the log of
#' the cdf
#' @export
cdf.series_system <- function(x) {

    R <- surv(x)
    p <- sum(nparams(x))

    function (t, theta = NULL, lower.limit = TRUE, log = FALSE) {

        theta <- params(x, theta)
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
#' @param object series_system object
#' @param log logical, log of survival function
#' @export
surv.series_system <- function(x) {

    R <- surv(x)
    p <- sum(nparams(x))

    function(t, theta = NULL, log = FALSE) {

        theta <- params(x, theta)
        stopifnot(length(theta) == p)

        R(t, theta, log = log)
    }
}

#' Hazard function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparams[j]} parameters.
#' @param h list of hazard functions
#' @export
hazard.series_system <- function(x) {

    h <- hazard(x)
    p <- sum(nparams(x))

    function(t, theta = NULL, log = FALSE) {
        theta <- params(x, theta)
        stopifnot(length(theta) == p)
        h(t, theta, log = log)
    }
}


#' cumulative hazard function for a component hazard function
#' @param haz hazard function
#' @export
cum_haz <- function(haz) {
    function(t, ...) {
        integrate(haz, lower = 0, upper = t, ...)$value
    }
}




















hazard_general_series_helper <- function(
    hazards,
    param_counts) {

    m <- length(hazards)
    function(t, theta, log.p = FALSE) {
        i0 <- 1
        i1 <- param_counts[1]
        v <- 0
        for (j in 1:m) {
            v <- v + hazards[[j]](t,
                theta[i0:i1],
                log.p = FALSE)
            i0 <- i1 + 1
            i1 <- i1 + param_counts[j]
        }
        ifelse(log.p, log(v), v)
    }
}

survival_general_series_helper <- function(
    survivals,
    param_counts) {

    m <- length(survivals)
    function(t, theta, log.p = FALSE) {
        i0 <- 1
        i1 <- param_counts[1]
        p <- ifelse(log.p, 0, 1)
        for (j in 1:m) {
            P <- survivals[[j]](t, theta[i0:i1], log.p = log.p)
            p <- ifelse(log.p, p + P, p * P)
            i0 <- i1 + 1
            i1 <- i1 + param_counts[j]
        }
        p
    }
}
