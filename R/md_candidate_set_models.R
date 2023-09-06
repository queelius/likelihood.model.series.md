#' Bernoulli candidate set model is a particular type of *uninformed* model.
#' Note that we do not generate candidate sets with this function. See
#' `md_cand_sampler` for that.
#'
#' This model satisfies conditions C1, C2, and C3.
#' The failed component will be in the corresponding candidate set with
#' probability 1, and the remaining components will be in the candidate set
#' with probability `p`.
#'
#' @param md masked data.
#' @param p a vector of probabilities
#' @param compvar column name of the component lifetime variables, defaults to
#'                `t`, e.g., `t1, t2, ..., tm`.
#' @param qvar column prefix for component probabilities, defaults to `q`,
#'             e.g., `q1, q2, ..., qm`.
#' @param deltavar column name of the right-censoring indicator variable, 
#'                 defaults to `delta`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @export
md_bernoulli_cand_c1_c2_c3 <- function(md, p,
    compvar = "t",
    qvar = "q",
    deltavar = "delta")
{
    n <- nrow(md)
    if (n == 0) {
        return(md)
    }
    p <- rep(p, length.out = n)
    Tm <- md_decode_matrix(md, compvar)
    if (is.null(Tm)) {
        stop("No component lifetime variables found")
    }
    m <- ncol(Tm)
    Q <- matrix(rep(p, m), nrow = n)
    Q[cbind(1:n, apply(Tm, 1, which.min))] <- 1
    if (!is.null(deltavar) && deltavar %in% colnames(md)) {
        Q[which(md[[deltavar]]),] <- 0
    }
    # remove in case it already has columns for q1,...,qm
    md[ ,paste0(qvar, 1:m)] <- NULL
    md %>% bind_cols(md_encode_matrix(Q, qvar)) %>%
           md_mark_latent(paste0(qvar, 1:m))
}

#' md_cand_sampler
#'
#' Candidate set generator. Requires columns for component probabilities
#' e.g., `q1,...,qm` where `qj` is the probability that the jth component
#' will be in the corresponding candidate set generated for that observation
#' in the `md` table.
#'
#' @param md masked data.
#' @param qvar column prefix for component probabilities, defaults to `q`,
#'             e.g., `q1, q2, ..., qm`.
#' @param setvar column prefix for candidate sets (as Boolean matrix), defaults
#'               to `x`, e.g., `x1, x2, ..., xm`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @importFrom stats runif
#' @export
md_cand_sampler <- function(md,qvar="q",setvar="x")
{
    Q <- md_decode_matrix(md,qvar)
    m <- ncol(Q)
    n <- nrow(Q)
    stopifnot(n > 0, m > 0)

    X <- matrix(rep(NA, m*n), nrow=n)
    for (i in 1:n)
        X[i,] <- runif(m) <= Q[i,]
    md %>% bind_cols(md_encode_matrix(X,setvar))
}

