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


series_system_from_names <- function(dist_names, param_counts) {
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

#' General series system object constructor.
#'
#' @param hazards list of hazard functions.
#' @param survivals list of reliability functions.
#' @param param_counts integer vector, we have m components, each of which
#'                     may have a different number of parameters. the j-th
#'                     component has `parma_counts[j]` parameters.
#' @return series_system object
#' @export
series_system <- function(
    hazards,
    survivals,
    samplers = NULL,
    param_counts) {
        stopifnot(length(hazards) == length(survivals),
                  length(hazards) == length(param_counts))

        # TODO: use `algebraic.dist` package
        structure(
            list(hazards = hazards,
                 survivals = survivals,
                 samplers = samplers,
                 param_counts = param_counts),
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
    stopifnot(is_series_system(object))
    h <- hazard_series_helper(object$hazards,
                              object$param_counts)
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
    stopifnot(is_series_system(object))
    S <- survival_series_helper(object$survivals,
                                object$param_counts)

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

    stopifnot(is_series_system(object))

    S <- survival(object)
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
    stopifnot(is_series_system(object))
    R <- survival(object)

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
    stopifnot(is_series_system(object))
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
    stopifnot(is_series_system(object))
    object$hazards
}

#' Number of parameters for a general series system.
#' @param object series_system object
#' @return number of parameters
#' @export
nparams.series_system <- function(object) {
    stopifnot(is_series_system(object))
    sum(object$param_counts)
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
qseries_system <- Vectorize(function(
    p,
    theta,
    object,
    options = list(), ...) {

        defaults <- list(
            alpha = 1, # initial step size, default is 1
            tol = 1e-3, # eps stopping condition, default is 1e-3
            t0 = 1,     # initial guess, default is 1
            ...)

        options <- modifyList(defaults, options)
        stopifnot(
            length(theta) == nparams(object),
            options$alpha > 0,
            options$tol > 0,
            options$t0 > 0,
            all(p > 0))

        h <- hazard(object)
        R <- survival(object)
        t1 <- NULL
        # TODO: replace with stats::optim
        repeat
        {
            r <- options$alpha
            repeat {
                t1 <- t0 - r * (1 - R(t, theta)) /
                    (h(t, theta) * R(t, theta))
                if (t1 > 0) break
                r <- r / 2
            }
            if (abs(t1 - t0) < options$tol) {
                break
            }
            t0 <- t1
        }
        t2
    }, vectorize.args = "p")

#' Sample from a general series system.
#'
#' @param n sample size
#' @param theta parameter vector, if NULL uses the values in `object`
#' @param object a `series_system` object
#' @param options list of options
#' @param ... additional arguments to pass into `options`
#' @export
rseries_system <- function(n, theta = NULL, object, options, ...) {

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

#' pdf for general series system
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param object series_system object
#' @export
dseries_system <- Vectorize(function(t, theta, object, log = FALSE) {
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    h <- hazard(object)
    R <- survival(object)
    if (is.null(theta))
        theta <- params(object)
    
    if (log)
        return(h(t,theta,log=log) - R(t,theta,log=log))
    else
        return(h(t,theta,log=log) / R(t,theta,log=log))
    
}, vectorize.args="t")

#' cdf for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param object series_system object
#' @export
pseries_system <- Vectorize(function(t, theta, object)
{
    R <- surv(object)
    1 - R(t, theta)
}, vectorize.args="t")

#' Survival function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param object series_system object
#' @param log logical, log of survival function
#' @export
surv_series_system <- Vectorize(function(
    t, theta = NULL, object, log = FALSE) {
    R <- surv(object)
    R(t, theta, log = log)
}, vectorize.args="t")

#' Hazard function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparams[j]} parameters.
#' @param h list of hazard functions
#' @export
hazard_series <- Vectorize(function(t, theta = NULL, object, log = FALSE) {
    h <- hazard(object)
    h(t,theta,log=log)
}, vectorize.args = "t")










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
