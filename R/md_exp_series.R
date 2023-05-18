#' Exponential series MLE
#' 
#' Maximum likelihood estimation functions for
#' Exponential series systems from masked data. Functions
#' include the log-likelihood, score, and FIM functions.
#' 
#' Masked component data approximately satisfies the
#' following conditions:
#' 
#' C1: Pr(K in C) = 1
#' C2: Pr(C=c | K=j, T=t) = Pr(C=c | K=j', T=t)
#'     for any j, j' in c.
#' C3: masking probabilities are independent of theta

#'
#' @author Alex Towell
#' @name Exponential series MLE
#' @keywords exponential, distribution, series, statistics, masked data
#' @seealso \code{\link{md_loglike_exp_series_C1_C2_C3}},
#'          \code{\link{md_score_exp_series_C1_C2_C3}},
#'          \code{\link{md_mle_exp_series_C1_C2_C3}}
#'          \code{\link{md_fim_exp_series_C1_C2_C3}}
NULL

#' Maximum likelihood estimator for exponential series
#' in which each component (or competing risk) is
#' parameterized by a single rate parameter.
#'
#' The sample (at least approximately) satisfy
#' competing risks conditions C1, C2, and C3 and
#' is failure time data may be right censored.
#'
#' @param md right-censored masked data
#' @param theta0 initial value for the MLE
#' @param control list of control parameters
#'  - eps convergence tolerance (default is 1e-6)
#'  - keep_obs logical, keep observations if TRUE
#'  - maxit integer, maximum number of iterations
#'  - compute_converged logical, whether to compute convergence info
#' @param ... pass additional arguments to `control`
#' @return MLE of type `md_mle_exp_series_C1_C2_C3`
#' 
#' @importFrom algebraic.mle mle_numerical
#' @importFrom algebraic.mle sim_anneal
#' @importFrom algebraic.mle newton_raphson
#' @importFrom MASS ginv
#' @importFrom algebraic.mle mle
#' @export
md_mle_exp_series_C1_C2_C3 <- function(md, theta0, control = list(), ...) {
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

    ll <- md_loglike_exp_series_C1_C2_C3(md, control$setvar, control$sysvar)
    if (control$use_simulated_annealing) {
        start <- sim_anneal(theta0, ll, control)
        theta0 <- start$par
    }

    sol <- NULL
    if (control$method == "newton-raphson") {
        sol <- newton_raphson(par = theta0, fn = ll, control)
    }
    else {
        sol <- optim(par = theta0, fn = ll, method = control$method, control)
    }

    if (sol$convergence != 0) {
        warning("optimization did not converge")
    }

    mle_sol <- mle_numerical(sol, control)
    class(mle_sol) <- c("md_mle_exp_series_C1_C2_C3", class(mle_sol))
    mle_sol
}

#' Generates a log-likelihood for an exponential series system with respect to
#' parameter `theta` for masked data with candidate sets that satisfy conditions
#' C1, C2, and C3.
#'
#' @param md masked data
#' @param setvar prefix of Boolean matrix encoding of candidate sets, defaults
#'               to `x`, e.g., `x1,...,xm`.
#' @param sysvar system lifetime (optionally right-censored) column name
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_loglike_exp_series_C1_C2_C3 <- function(md,setvar="x",sysvar="t")
{
    stopifnot(sysvar %in% colnames(md))
    sum.t <- sum(md[,sysvar])
    if ("delta" %in% colnames(md))
        md <- md %>% filter(.data$delta == FALSE)
    md$C <- md_decode_matrix(md, setvar)
    md <- md %>% group_by(.data$C) %>% count()
    n <- nrow(md)

    function(theta)
    {
        if (any(theta <= 0)) return(NA)
        f <- -sum.t * sum(theta)
        for (i in 1:n)
            f <- f + md$n[i] * log(sum(theta[md$C[i, ]]))
        f
    }
}


#' Generates a score function for an exponential series system (or competing risks)
#' with respect to parameter `theta` for masked component failure (or masked
#' competing risks) with candidate sets (risks) that approximately satisfy conditions
#' C1, C2, and C3.
#'
#' @param md right-censored failure-time data with masked competing risks
#' @param setvar prefix of Boolean matrix encoding of candidate sets, defaults
#'               to `x`, e.g., `x1,...,xm`.
#' @param sysvar system lifetime (optionally right-censored) column name
#' @return score function of type `R^m -> R`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_score_exp_series_C1_C2_C3 <- function(md,setvar="x",sysvar="t")
{
    stopifnot(sysvar %in% colnames(md))
    sum.t <- sum(md[,sysvar])
    if ("delta" %in% colnames(md))
        md <- md %>% filter(.data$delta == FALSE)
    md$C <- md_decode_matrix(md, setvar)
    md <- md %>% group_by(.data$C) %>% count()
    m <- ncol(md$C)

    function(theta)
    {
        if (any(theta <= 0)) return(NA)
        v <- rep(-sum.t, m)
        for (j in 1:m)
        {
            for (i in 1:nrow(md))
            {
                if (md$C[i, j]) {
                    v[j] <- v[j] + md$n[i] / sum(theta[md$C[i, ]])
                }
            }
        }
        as.matrix(v)
    }

    #function(theta)
    #{
    #    if (any(theta <= 0))
    #        return(NA)
    #    # Use sapply and vectorization to replace the nested for-loops
    #    v <- sapply(1:m, function(j) sum(md$n[md$C[, j]] / rowSums(theta[md$C])))
    #    as.matrix(v - sum.t)
    #}    
}

#' Generates the observed information matrix (FIM) for an exponential series system
#' with respect to parameter `theta` for masked data with candidate sets
#' that approximately satisfy conditions C1, C2, and C3.
#'
#' @param md masked data with candidate sets that meet the
#'           regular candidate model
#' @param setvar prefix of Boolean matrix encoding of candidate sets, defaults
#'               to `x`, e.g., `x1,...,xm`.
#' @param sysvar system lifetime (optionally right-censored) column name
#' @return observed information matrix of type `R^m -> R^(m x m)`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_fim_exp_series_C1_C2_C3 <- function(md, setvar="x", sysvar="t")
{
    if ("delta" %in% colnames(md))
        md <- md %>% filter(.data$delta == FALSE)
    md$C <- md_decode_matrix(md, setvar)
    md <- md %>% group_by(.data$C) %>% count()
    m <- ncol(md$C)
    n <- nrow(md)

    function(theta)
    {
        if (any(theta <= 0)) return(NA)
        I <- matrix(rep(0, m * m), nrow = m)
        for (j in 1:m)
        {
            for (k in 1:m)
            {
                for (i in 1:nrow(md))
                {
                    if (md$C[i, j] && md$C[i, k]) {
                        I[j, k] <- I[j, k] + md$n[i] / sum(theta[md$C[i, ]])^2
                    }
                }
            }
        }
        I
    }

    #function(theta)
    #{
    #    if (any(theta <= 0))
    #        return(limit)
    #    # Use outer and sapply to replace the nested for-loops
    #    I <- outer(1:m, 1:m, function(j, k) {
    #        ifelse(md$C[, j] & md$C[, k], 
    #               sum(md$n / rowSums(theta[md$C])^2), 
    #               0)
    #    })
    #    I
    #}    
}
