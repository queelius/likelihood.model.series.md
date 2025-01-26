#' Bernoulli candidate set model for systems with unobserved components.
#' 
#' Bernoulli candidate set model is a particular type of *uninformed* model.
#' Note that we do not generate candidate sets with this function. See
#' `md_cand_sampler` for that.
#'
#' This model satisfies conditions C1, C2, and C3.
#' The failed component will be in the corresponding candidate set with
#' probability 1, and the remaining components will be in the candidate set
#' with probability `p` (the same probability for each component). `p` 
#' may be different for each system, but it is assumed to be the same for
#' each component within a system, so `p` can be a vector such that the
#' length of `p` is the number of systems in the data set (with recycling
#' if necessary).
#'
#' @param df masked data.
#' @param p a vector of probabilities (p[j] is the probability that the jth
#'          system will include a non-failed component in its candidate set,
#'          assuming the jth system is not right-censored).
#' @param prob column prefix for component probabilities, defaults to
#'             `q`, e.g., `q1, q2, q3`.
#' @param comp column prefix for component lifetimes, defaults to `t`,
#'             e.g., `t1, t2, t3`.
#' @param right_censoring_indicator right-censoring indicator column name.
#'     if TRUE, then the system lifetime is right-censored, otherwise it is
#'     observed. If NULL, then no right-censoring is assumed. Defaults to
#'     `delta`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix md_mark_latent
#' @importFrom dplyr %>% bind_cols
#' @export
md_bernoulli_cand_c1_c2_c3 <- function(
    df,
    p,
    prob = "q",
    comp = "t",
    right_censoring_indicator = "delta") {

    n <- nrow(df)
    if (n == 0) {
        return(df)
    }
    p <- rep(p, length.out = n)
    Tm <- md_decode_matrix(df, comp)
    if (is.null(Tm)) {
        stop("No component lifetime variables found")
    }
    m <- ncol(Tm)
    Q <- matrix(p, nrow = n, ncol = m)
    Q[cbind(1:n, apply(Tm, 1, which.min))] <- 1

    if (!is.null(right_censoring_indicator)) {

        if (!right_censoring_indicator %in% colnames(df)) {
            stop("right_censoring_indicator not in colnames(df)")
        }
        Q[!df[[right_censoring_indicator]], ] <- 0
    }

    # remove in case it already has columns for q1,...,qm
    df[ , paste0(prob, 1:m)] <- NULL
    df %>% bind_cols(md_encode_matrix(Q, prob)) %>%
           md_mark_latent(paste0(prob, 1:m))
}

#' Sample candidate sets for systems with unobserved components.
#' 
#' Candidate set generator. Requires columns for component probabilities
#' e.g., `q1,...,qm` where `qj` is the probability that the jth component
#' will be in the corresponding candidate set generated for that observation
#' in the `md` table.
#'
#' @param df (masked) data frame
#' @param prob column prefix for component probabilities, defaults to
#'             `q`, e.g., `q1, q2, q3`.
#' @param candset column prefix for candidate sets (as Boolean matrix),
#'                defaults to `x`, e.g., `x1, x2, x3`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @importFrom stats runif
#' @export
md_cand_sampler <- function(df, prob = "q", candset = "x") {

    if (!is.data.frame(df)) {
        stop("df must be a data frame")
    }

    n <- nrow(df)
    if (n == 0) {
        return(df)
    }

    Q <- md_decode_matrix(df, prob)
    m <- ncol(Q)
    if (m == 0) {
        stop("No component probabilities found")
    }

    X <- matrix(NA, nrow=n, ncol=m)
    for (i in 1:n) {
        X[i, ] <- runif(m) <= Q[i, ]
    }
    df %>% bind_cols(md_encode_matrix(X, candset))
}

