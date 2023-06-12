loglik.exp_right_censored = function(row, rates) {
    -sum(rates) * row$t
}

loglik.exp_exact_masked_component = function(row, rates) {
  X <- select(row, starts_with("x"))
  -sum(rates) * row$t + log(sum(rates[X]))
}

score.exp_right_censored = function(row, rates) {
  rep(-row$t, length(rates))
}

score.exp_exact_masked_component = function(row, rates) {
    X <- select(row, starts_with("x"))
    X / sum(rates[X]) - row$t
}

hess_loglik.exp_right_censored = function(row, rates) {
  m <- length(rates)
  matrix(rep(0, m * m), nrow = m, ncol = m)
}
  
exact_failure_time_masked_component = function(row, rates) {
  X <- select(row, starts_with("x"))
  outer(X, X, function(x, y) ifelse(x & y, 1 / sum(rates[X])^2, 0))
}






###########



#' Likelihood contribution for an exponential series system with respect to
#' rate parameters for masked component cause of failure with candidate sets
#' that satisfy conditions C1, C2, and C3.
#'
#' @param t system failure time
#' @param X candidate set (Boolean vector)
#' @param rates rate parameter
#' @return likelihood contribution
#' @export
loglik.exp_series_md_C1_C2_C3 <- function(t, X, rates, log = TRUE) {
    stopifnot(length(rates) > 0, length(rates) == length(X), all(rates > 0))
    if (log) {
        -t * sum(theta) + log(sum(rates[X]))
    } else {
        exp(-t * sum(theta)) * sum(rates[X])
    }
}

#' Gradient of the log-likelihood contribution for an exponential series system 
#' with respect to rate parameters for masked component cause of failure with 
#' candidate sets that satisfy conditions C1, C2, and C3.
#'
#' @param t system failure time
#' @param X candidate set (Boolean vector)
#' @param rates rate parameters
#' @return score of likelihood contribution
#' @export
score_exp_series_md_C1_C2_C3 <- function(t, X, rates) {
    stopifnot(length(rates) > 0, length(rates) == length(X), all(rates > 0))
    -t + X / sum(rates[X])
}

#' Hessian of the log-likelihood contribution (negative of the
#' observed FIM) for an exponential series system with respect to rate
#' parameters for masked component cause of failure with candidate sets
#' that satisfy conditions C1, C2, and C3.
#'
#' @param t system failure time
#' @param X candidate set (Boolean vector)
#' @param rates rate parameters
#' @return hessian of likelihood contribution
#' @export
hess_exp_series_md_C1_C2_C3 <- function(t, X, rates) {
    m <- length(rates)
    stopifnot(m > 0, m == length(X), all(rates > 0))
    J <- matrix(rep(0, m * m), nrow = m)
    for (j in 1:m) {
        for (k in 1:j) {
            if (X[j] && X[k]) {
                J[j, k] <- 1 / sum(rates[X])^2
                J[k, j] <- J[j, k]
            }
        }
    }
    J
}














##################
#' Fisher information matrix for exponential series system under
#' a candidate set model that satisfies conditions C1, C2, and C3
#' called the Bernoulli candidate set model.
#' 
#' MC simulation of the FIM for a single observation
#' of the exponential series system with masked
#' component failure that follows the Bernoulli candidate set
#' model (which satisfies conditions C1, C2, and C3) where:
#' 
#'   (1) for a given observation, `p ~ rp`, and `p` is the
#'       probability that each of the non-failed components
#'       are in the corresponding candidate set.
#' 
#'   (2) for a given observation, `tau ~ rtau`, and `tau` is
#'       the observable right-censoring time.
#' 
#' The data is generated from the true `theta` and the FIM is
#' estimated by averaging over `R` simulations.
#' 
#' @param rates vector of component failure rates
#' @param rtau function that generates right-censoring times
#'            for a single observation
#' @param rp function that generates `p` for a single observation
#' @param R number of simulations
#' @export
md_expected_fim_exp_series_system_C1_C2_C3 <- function(
    rates,
    rtau,
    rp,
    R = 999)
{
    m <- length(rates)
    if (m == 0) {
        stop("Must have at least one component")
    }
    if (any(rates <= 0)) {
        stop("`rates` must be positive")
    }

    fim <- matrix(0, nrow = m, ncol = m)
    tau <- rtau(R)
    p <- rp(R)
    for (i in 1:R) {
        data <- md_exp_series_system_bernoulli_cand_C1_C2_C3(
            n = 1,
            rates = rates,
            p = p[i],
            tau = tau[i])

        v <- md_score_exp_series_C1_C2_C3(data)(rates)
        fim <- fim + v %*% t(v)
    }
    fim / R
}

exp_fim_2 <- function(rates, rtau, rp, R = 999) {
    m <- length(rates)
    if (m == 0) {
        stop("Must have at least one component")
    }
    if (any(rates <= 0)) {
        stop("`rates` must be positive")
    }

    fim <- matrix(0, nrow = m, ncol = m)
    tau <- rtau(R)
    p <- rp(R)
    data <- exp_series_md_bernoulli_cand_C1_C2_C3(
        n = R,
        rates = rates,
        p = p,
        tau = tau)

    for (i in 1:R) {
        v <- md_score_exp_series_C1_C2_C3(data[i,])(rates)
        fim <- fim + v %*% t(v)
    }
    fim / R
}

exp_fim_3 <- function(rates, rtau, rp, P, R = 999) {
    m <- length(rates)
    if (m == 0) {
        stop("Must have at least one component")
    }
    if (any(rates <= 0)) {
        stop("`rates` must be positive")
    }

    fim <- matrix(0, nrow = m, ncol = m)
    tau <- rtau(R)
    p <- rp(R)
    data <- exp_series_md_bernoulli_cand_C1_C2_C3(
        n = R,
        rates = rates,
        p = p,
        tau = tau)

    for (i in 1:R) {
        v <- md_score_exp_series_C1_C2_C3(data[i,])(rates)
        fim <- fim + v %*% t(v)
    }
    fim / R
}




#' Hessian of the log-likelihood contribution (negative of the
#' observed FIM) for an exponential series system with respect to rate
#' parameters for masked component cause of failure with candidate sets
#' that satisfy conditions C1, C2, and C3.
#'
#' @param t system failure time
#' @param X candidate set (Boolean vector)
#' @param rates rate parameters
#' @return hessian of likelihood contribution
#' @export
hess_exp_series_md_C1_C2_C3_2 <- function(t, X, rates) {
    m <- length(rates)
    stopifnot(m > 0, m == length(X), all(rates > 0))
    outer(X, X, function(x, y) ifelse(x & y, 1 / sum(rates[X])^2, 0))
}
