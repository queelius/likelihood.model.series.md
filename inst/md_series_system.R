#' Likelihood model for series systems defined by hazard functions and
#' survival functions, for masked component cause of failure (candidate
#' sets) satisfying conditions C1, C2, and C3, and right censoring where
#' the censoring indicator is independent of the cause of failure and
#' the time of failure.
#'
#' @author Alex Towell
#' @name General series likelihood model for masked data
#' @keywords distribution, series, statistics, masked data
#' @seealso \code{\link{md_loglike_general_series_C1_C2_C3}},
#'          \code{\link{md_score_general_series_C1_C2_C3}},
#'          \code{\link{md_mle_general_series_C1_C2_C3}}

likelihood_md_series_system_C1_C2_C3 <- function(x) {
    structure(list(x = x),
        class = c("likelihood_md_series_system_C1_C2_C3",
                  "likelihood_model"))
}

#' Masked data approximately satisfies the following set of conditions:
#' C1: Pr(K[i] in C[i]) = 1
#' C2: Pr(C[i]=c[i] | K[i]=j, T[i]=t[i]) =
#'         Pr(C[i]=c[i] | K[i]=j', T[i]=t[i])
#'     for any j,j' in c[i].
#' C3: masking probabilities are independent of theta
#'
#' @param md masked data
#' @param nparams integer vector, we have m components, each of which may have a
#'                different number of parameters. the j-th component has
#'                \code{nparmas[j]} parameters.
#' @param h list of hazard functions
#' @param R list of survival functions, defaults to NULL, in which case we
#'          generate the survival functions from \code{h}.
#' @returns a log-likelihood function with respect to theta given \code{md}
#' @importFrom md.tools md_decode_matrix
#' @export
loglik.likelihood_md_series_system_C1_C2_C3 <- function(x) {
    
    
    function(md, theta) {

    m <- length(x$hazards)

    if (is.null(R)) {
        R <- list(length=m)
        j <- 1
        for (hj in h) {
            R[[j]] <- generate_survival_from_hazard(hj)
            j <- j + 1
        }
    }
    stopifnot(m > 0, m == length(R))
    #series_h <- hazard_general_series_helper(h,nparams)
    #series_R <- survival_general_series_helper(R,nparams)

    C <- md_decode_matrix(md,"x")
    stopifnot(ncol(C) == m)
    n <- nrow(md)
    stopifnot(n > 0)

    right_censoring <- "delta" %in% colnames(md)
    stopifnot(!right_censoring || "s" %in% colnames(md))
    t <- ifelse(right_censoring, md$s, md$t)

    log.R <- function(t, theta) {
        sum <- 0
        for (Rj in R) {
            sum <- sum + log(Rj(t, theta))
        }
        sum
    }

    # total number of params
    p <- sum(nparams)

    function(theta)
    {
        stopifnot(length(theta) == p)
        res <- 0
        for (i in 1:n) {
            if (right_censoring && md$delta[i]) {
                res <- res + log.R(t[i], theta)
            } else { #if (!right_censoring || (right_censoring && !md$delta[i]))
                haz <- 0
                for (j in (1:m)[C[i,]])
                    haz <- haz + h[[j]](t[i],theta)
                res <- res + log(haz) * log.R(t[i], theta)
            }
        }
        res
    }
}
