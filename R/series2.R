#' General series
#'
#' This file contains functions related to a general series distribution.
#' Functions include simulation, pdf, cdf, quantile, and other related
#' functions.
#' 
#' For parameter estimation from masked data, see the file \code{md_series_mle.R}.
#'
#' @author Alex Towell
#' @name General series
#' @keywords distribution, series, statistics
NULL



#'
#' @param hazards list of hazard functions, if NULL, then survivals must be
#'                provided and hazard functions will be generated from them
#' @param param_counts integer vector, we have m components, each of which
#'                     may have a different number of parameters. the j-th
#'                     component has `parma_counts[j]` parameters.
#' @param survivals list of reliability functions, default is NULL (numerical
#'                  approximation from hazard functions will be used)
#' @param options list of options
series_system <- function(
    hazards,
    param_counts,
    survivals = NULL,
    options = list(), ...) {

        defaults <- list(
            pdfs = NULL,
            samplers = NULL,
            cdfs = NULL,
            quantiles = NULL,
            ...)

        options <- modifyList(defaults, options)

        m <- length(hazards)
        stopifnot(m == length(param_counts))

        if (is.null(survivals)) {
            survivals <- list()
            for (i in 1:m)
                survivals[[i]] <- generate_survival_from_hazard(
                    hazards[[i]])
        }
        stopifnot(m == length(survivals))

        structure(
            list(
                hazards = hazards,
                survivals = survivals,
                pdfs = options$pdfs,
                samplers = options$samplers,
                cdfs = options$cdfs,
                quantiles = options$quantiles,
                param_counts = param_counts),
            class = c("series_system", "distribution"))
}


#' @export
hazards.series_system <- function(object, ...) {
    object$hazards
}

#' @export
nparams.series_system <- function(object, ...) {
    sum(object$param_counts)
}



#' Quantile function (inverse of the cdf) for a general series system.
#' By definition, the quantile \code{p * 100%} is the value \code{t} that
#' satisfies \code{F(t) - p = 0}. We solve for \code{t} using Newton's method.
#'
#' @param p vector of probabilities.
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of reliability functions, default is NULL (numerical approximation from h)
#' @export
qseries_system <- Vectorize(function(
    p,
    theta,
    object,
    options = list(), ...) {

        defaults <- list(
            alpha0 = 1, # initial step size, default is 1
            eps = 1e-3, # eps stopping condition, default is 1e-3
            t0 = 1)     # initial guess, default is 1

        options <- modifyList(defaults, options)
        stopifnot(length(theta) == nparams(object),
            alpha0 > 0,
            option$eps > 0,
            all(p > 0))

        h <- hazard_series_helper(object$hazards, object$param_counts)
        R <- survival_series_helper(object$survivals, object$param_counts)
        t1 <- NULL
        repeat
        {
            r <- options$alpha0
            repeat {
                t1 <- t0 - r * (1 - R(t, theta)) /
                    (h(t, theta) * R(t, theta))
                if (t1 > 0) break
                r <- r / 2
            }
            if (abs(t2 - t1) < options$eps) {
                break
            }
            t1 <- t2
        }
        t2
    }, vectorize.args = "p")

#' Sample from a general series system of \code{m} components whose hazard
#' and reliability functions are respectively given by \code{h} and \code{R}.
#'
#' @param n sample size
#' @param theta parameter vector
#' @param object a `series` system object
#' @export
rseries <- function(n, theta, object) {

    stopifnot(is_series_object(object))

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
    
    system_lifetime <- apply(lifetimes, 1, min)
    if (include_latent) {
        md.tools::md_encode_matrix(lifetimes, "t") %>% 
            md.tools::md_encode_matrix(system_lifetime, "t") %>% 
            md.tools::md_encode_matrix(comps, "t") %>% 
            md.tools::md_encode_matrix(theta, "t") %>% 
            md.tools::md_encode_matrix(nparams, "t") %>% 
            md.tools::md_encode_matrix(h, "t") %>% 
            md.tools::md_encode_matrix(R, "t") %>% 
            md.tools::md_encode_matrix(include_latent, "t") -> latent
    } else {
        return(system_lifetime)
    }    
}

#' pdf for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of reliability functions
#' @export
dseries <- Vectorize(function(t,theta,nparams,h,R)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    series_h <- hazard_series_helper(h,nparams)
    series_R <- survival_series_helper(R,nparams)
    series_h(t,theta) / series_R(t,theta)
}, vectorize.args="t")

#' cdf for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param R list of reliability functions
#' @export
pseries <- Vectorize(function(t,theta,nparams,R)
{
    stopifnot(length(theta)==sum(nparams))
    stopifnot(length(h)==length(R))
    series_R <- survival_series_helper(R,nparams)
    1-series_R(t,theta)
}, vectorize.args="t")

#' Survival function for general series
#'
#' @param t series system lifetime
#' @param theta parameter vector
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparams[j]} parameters.
#' @param R list of reliability functions
#' @note convert this into a generator over (t,theta)
#' @export
survival_series <- Vectorize(function(t,theta,nparams,R)
{
    stopifnot(length(theta)==sum(nparams))
    series_R <- survival_series_helper(R,nparams)
    function(t,theta) series_R(t,theta)
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
hazard_series <- Vectorize(function(t,theta,nparams,h)
{
    stopifnot(length(theta)==sum(nparams))
    series_h <- hazard_series_helper(h,nparams)
    series_h(t,theta)
}, vectorize.args="t")

