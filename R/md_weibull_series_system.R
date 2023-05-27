#' Maximum likelihood estimation functions for Weibull series systems
#' from masked data.
#'
#' Functions include the log-likelihood, score, and MLE functions.
#' 
#' Masked component data approximately satisfies the following conditions:
#' C1: Pr(K in C) = 1
#' C2: Pr(C=c | K=j, T=t) = Pr(C=c | K=j', T=t)
#'     for any j, j' in c.
#' C3: masking probabilities are independent of theta
#'
#' @author Alex Towell
#' @name Weibull series MLE
#' @keywords weibull, distribution, series, statistics, masked data
#' @seealso \code{\link{md_loglike_weibull_series_C1_C2_C3}},
#'          \code{\link{md_score_weibull_series_C1_C2_C3}},
#'          \code{\link{md_mle_weibull_series_C1_C2_C3}}
NULL


#' Maximum likelihood estimator for Weibull series
#' in which each component (or competing risk) is
#' parameterized by a shape and scale parameter (k and lambda).
#'
#' The sample (at least approximately) satisfy
#' competing risks conditions C1, C2, and C3 and
#' its system failure time data may be right censored.
#'
#' @param md right-censored masked data
#' @param theta0 initial value for the MLE
#' @param control list of control parameters
#'  - eps convergence tolerance (default is 1e-6)
#'  - keep_obs logical, keep observations if TRUE
#'  - maxit integer, maximum number of iterations
#'  - compute_converged logical, whether to compute convergence info
#' @param ... pass additional arguments to `control`
#' @return MLE of type `md_mle_weibull_series_C1_C2_C3`
#' @importFrom algebraic.mle mle_numerical sim_anneal mle
#' @importFrom MASS ginv
#' @export
md_mle_weibull_series_C1_C2_C3 <- function(md, theta0, control = list(), ...) {
    defaults <- list(
        keep_obs = FALSE,
        maxit = 1000L,
        compute_converged = FALSE,
        use_simulated_annealing = TRUE,
        fnscale = -1,
        sysvar = "t",
        setvar = "x",
        method = "BFGS",
        lower_bound = 1e-6, # should all be positive
        sup = function(theta) all(theta >= control$lower_bound),
        proj = function(theta) pmax(theta, control$lower_bound))

    control <- modifyList(defaults, control)
    control <- modifyList(control, list(...))

    ll <- md_loglike_weibull_series_C1_C2_C3(md, control$setvar, control$sysvar)
    if (control$use_simulated_annealing) {
        start <- sim_anneal(theta0, ll, control)
        theta0 <- start$par
    }
    sol <- optim(
        par = theta0,
        fn = ll, hessian = TRUE,
        method = control$method,
        control)

    if (sol$convergence != 0) {
        warning("optimization did not converge")
    }

    mle_numerical(
        sol = sol,
        options = control,
        superclasses = c(
            "md_mle_weibull_series_C1_C2_C3",
            "md_mle_weibull_series",
            "md_mle_series",
            "md_mle"))
}

#' Log-likelihood generator for Weibull series system on masked data,
#' where the masked data is in the form of right-censored system lifetimes
#' and masked component cause of failure.
#'
#' @param md masked data
#' @returns a log-likelihood function with respect to theta given md
#' @importFrom md.tools md_decode_matrix
#' @export
md_loglike_weibull_series_C1_C2_C3 <- function(md)
{
    right_censoring <- "delta" %in% colnames(md)
    t <- NULL
    delta <- NULL
    if (right_censoring)
    {
        stopifnot("s" %in% colnames(md))
        t <- md$s
        delta <- md$delta
    }
    else
    {
        stopifnot("t" %in% colnames(md))
        t <- md$t
    }

    C <- md_decode_matrix(md,"x")
    m <- ncol(C)
    n <- nrow(md)
    stopifnot(m > 0)
    stopifnot(n > 0)
    md <- NULL

    function(theta)
    {
        # theta should be a parameter vector of length 2*m
        scales <- theta[(0:(m-1)*2)+1]
        shapes <- theta[(1:m)*2]

        s <- 0
        for (i in 1:n)
        {
            s <- s - sum((t[i]/scales)^shapes)
            if (!right_censoring || !delta[i])
                s <- s + log(sum(shapes[C[i,]]/scales[C[i,]]*(t[i]/scales[C[i,]])^(shapes[C[i,]]-1)))
        }
        s
    }
}

#' Score function generator (gradient of the log-likelihood) for Weibull series
#' system on masked data, where the masked data is in the form of right-censored
#' system lifetimes and masked component cause of failure. Note that if the
#' right censoring indicator variable `delta` is missing, we assume that all
#' system lifetimes are observed.
#'
#' Masked component data `md` approximately satisfies the following conditions:
#' C1: Pr(K in C) = 1
#' C2: Pr(C=c | K=j, T=t) = Pr(C=c | K=j', T=t)
#'     for any j, j' in c.
#' C3: masking probabilities are independent of theta
#'
#' @param md masked data
#' @param sysvar name of the system lifetime column
#' @param setvar prefix symbol for the candidate sets column
#' @returns a score function with respect to theta given md
#' @importFrom md.tools md_decode_matrix
#' @export
md_score_weibull_series_C1_C2_C3 <- function(md, sysvar="t", setvar="x")
{
    stopifnot(sysvar %in% colnames(md))
    C <- md_decode_matrix(md,setvar)
    stopifnot(!is.null(C))
    m <- ncol(C)
    n <- nrow(C)
    stopifnot(m > 0, n > 0)
    delta <- ifelse("delta" %in% colnames(md), md$delta, rep(FALSE,nrow(md)))
    t <- md[,sysvar] # system lifetimes
    md <- NULL

    function(theta)
    {
        # theta should be a parameter vector of length 2*m
        scales <- theta[(0:(m-1)*2)+1]
        shapes <- theta[(1:m)*2]

        scr <- rep(0,m)
        for (j in 1:m)
        {
            s <- 0
            for (i in 1:n)
            {
                s <- s - shapes[j]/scales[j]*(t[i]/scales[j])^shapes[j]/(sum((t[i]/scales)^shapes))
                if (C[i,j] && !delta[i])
                    s <- s - shapes[j]^2 * (t[i]/scales[j])^(shapes[j]-1)/
                        sum(shapes[C[i,]]/scales[C[i,]]*(t[i]/scales[C[i,]])^(shapes[C[i,]]-1))
            }
            scr[j] <- s
        }
    }
}

