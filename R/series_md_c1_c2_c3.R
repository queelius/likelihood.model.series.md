#' Likelihood model for series systems with masked component cause
#' of failure with candidate sets that satisfy conditions C1, C2, and C3,
#' desribed below. It is defined by hazard functions and
#' survival functions
#' 
#' This model satisfies the concept of a `likelihood_model` in the
#' `likelihood_model` package by providing the following methods:
#' 
#'  (1) `loglik.exp_series_md_c1_c2_c3`
#'  (2) `score.exp_series_md_c1_c2_c3`
#'  (3) `hess_loglik.exp_system_series_md_c1_c2_c3`
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
#' @author Alex Towell
#' @keywords hazard, survival, distribution, series, statistics, masked data
series_md_c1_c2_c3 <- function(x) {
    structure(list(x = x),
        class = c("series_md_c1_c2_c3",
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
#' @returns a log-likelihood function with respect to theta
#' @importFrom md.tools md_decode_matrix
#' @importFrom likelihood.model loglik
#' @export
loglik.series_md_c1_c2_c3 <- function(model) {

    function(df, theta, candset = "x", set.var = "t") {
        C <- md_decode_matrix(md, candset)
        if (is.null(C)) {
            stop("No candidate set variables found")
        }

        n <- nrow(df)
        stopifnot(n > 0)

        # log.R <- function(t, theta) {
        #     sum <- 0
        #     for (Rj in R) {
        #         sum <- sum + log(Rj(t, theta))
        #     }
        #     sum
        # }

        # total number of params
        p <- sum(x$nparams)
        stopifnot(length(theta) == p)
        res <- 0
        for (i in 1:n) {
            haz <- 0
            for (j in (1:m)[C[i,]]) {
                haz <- haz + h[[j]](t[i],theta)
            }
            res <- res + log(haz) * log.R(t[i], theta)
        }
        res
    }
}