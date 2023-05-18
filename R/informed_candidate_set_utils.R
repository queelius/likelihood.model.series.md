#' kl_divergence_bernoulli
#'
#' KL divergence between two Bernoulli random vectors.
#' p = (p1,p2,...,pn) and q = (q1,q2,...,qn).
#'
#' @param p probability of success for the first distribution
#' @param q probability of success for the second distribution
#' @return KL divergence between the two distributions
#' @export
#' @examples
#' kl_divergence_bernoulli(c(0.1,0.1), c(0.1,0.1)) = 0
kl_divergence_bernoulli <- function(p, q)
{
    stopifnot(length(p) == length(q))
    stopifnot(all(p >= 0 & p <= 1))
    stopifnot(all(q >= 0 & q <= 1))

    kl <- function(x, y)
    {
        if (x == 0 || x == 1 || y == 0 || y == 1)
        {
            if (x == y)
                return(0)
            else if (x == 0)
                return((1 - x) * log2((1 - x) / (1 - y)))
            else if (x == 1)
                return(x * log2(x / y))
        }
        return(x * log2(x / y) + (1 - x) * log2((1 - x) / (1 - y)))
    }
    sum(mapply(kl, p, q))
}

#' informative_masking_by_rank
#'
#' Returns `m=length(ts)` probabilities, where the j-th probability is the
#' probability that the j-th component will be in the candidate set.
#' probabilities are a function of rank (rather than by component failure
#' times).
#'
#' The shape of the masking is given by two parametrs, `alpha`, which must be
#' non-negative, and `beta`, which must be between 0 and 1.
#' As `alpha` goes to 0, the probability vector approaches
#'     (beta, beta, ..., 1, beta, ..., beta)
#' where the 1 is for the failed component and as `alpha` goes to infinity, it
#' assigns probability 1 and `beta` respectively to rank 1 (failed component)
#' and rank 2 components (and the remaining probabilities are zero).
#'
#' @param ts component failure times for the series system
#' @param alpha rate of change with respect to rank, as alpha -> 0 probabilities
#'              go to discrete uniform distribution
#' @param beta max weight (for ranked 2 component)
#' @export
informative_masking_by_rank <- function(ts, alpha, beta)
{
    stopifnot(alpha0 >= 0)
    stopifnot(beta0 >= 0 && beta0 <= 1)

    m <- length(ts)
    stopifnot(m > 0)

    ranks <- 0:(m-2)
    wts <- beta * exp(-alpha * ranks)
    probs <- c(1, wts)

    probs[order(ts)] <- probs
    probs
}

#' generate_relaxed_cand_C1
#'
#' Generates a probability `Q = (q1, q2, ..., qm)` such that `qj` is the
#' probability that the j-th component is in the candidate set, `qk = 1`, where
#' `k` is failed component. `Q` is an *informed* candidate model that uses
#' `informative_masking_by_rank` to assign higher probabilities to components
#' that failed earlier (which is something we typically only know in, say, a
#' simulation study).
#'
#' The probabilities `Q` have two constraints on them. Let
#' `P = (p, ..., p, 1, p, ..., p)` be the bernoulli candidate model that
#' satisfies conditions C1, C2, and C3. Then, the KL-divergence between `P`
#' and `Q` is as close as possible to `d` while satisfying `sum(P) == sum(Q)`.
#'
#' For `d = 0`, `Q == P`. As `d` increases, `Q` becomes more informative about
#' the components. Given the structure of `informative_masking_by_rank`, it may
#' not be possible to satisfy every `d` specified, but we get as close as
#' we can, which should permit useful experiments.
#'
#' @param d numeric, the KL divergence from P = (p, p, ..., p, 1, p, ..., p)
#'          to try to obtain
#' @param ts component failure times for the series system
#' @param p numeric, defines `P = (p, ..., p, 1, p, ..., p)`.
#' @param debug Logical, whether to output debugging information while running
#' @param eps numeric, stopping condition.
#' @param alpha0 numeric, initial guess for `alpha` parameter of
#'               `informative_masking_by_rank`.
#' @param beta0 numeric, initial guess for `beta` parameter of
#'              `informative_masking_by_rank`.
#' @param lambda numeric, controls how much the two constraints are weighted.
#'               Lower value specifies more enforcement of the KL-divergence
#'               constraint being closer to `d`. Defaults to 1.
#' @param max_iter Integer, maximum number of iterations before giving up.
#' @param lr numeric, learning rate.
#' @export
generate_relaxed_cand_C1 <- function(d, ts, p,
                                     debug=F,
                                     eps=1e-8,
                                     alpha0=1,
                                     beta0=p,
                                     lambda=1,
                                     max_iter=10000L,
                                     lr=1)
{
    stopifnot(alpha0 >= 0)
    stopifnot(beta0 >= 0 && beta0 <= 1)

    m <- length(ts)
    P = rep(p,m)
    k <- which.min(ts)
    P[k] <- 1
    if (d == 0)
        return(list(P=P,Q=P,KL=0,alpha=Inf,beta=p,converged=TRUE,ts=ts,k=k))

    P.sum <- sum(P)
    cost <- function(param) {
        alpha <- param[1]
        beta <- param[2]
        Q <- informative_masking_by_rank(ts, alpha, beta)
        delta <- sum(Q)-P.sum
        KL <- kl_divergence_bernoulli(P,Q)

        if (debug)
        {
            cat("-----------------------\n")
            cat("(alpha,,beta) = (", alpha, ",", beta, ")\n")
            cat("P = ", P, "\n")
            cat("Q = ", Q, "\n")
            #cat("P = ", P, "\n")
            cat("constraint 1: (sum(Q) - sum(P)) = ", delta^2, "\n")
            cat("KL = ", KL, "\n")
            cat("constraint 2: (KL - d)^2 = ", (KL - d)^2, "\n")
            cat("=======================\n")
        }

        (KL - d)^2 + lambda*delta^2
    }

    sup <- function(theta) theta[1] >= 0 && theta[2] >= 0 && theta[2] <= 1
    res <- grad_descent(f=cost, x0=c(alpha0,beta0),
                        sup=sup, debug=debug, lr=lr, eps=eps, max_iter=max_iter)
    if (!res$converged)
        return(list(Q=NULL,P=P,KL=NULL,alpha=NULL,beta=NULL,converged=FALSE,ts=ts,k=k))

    alpha <- res$param[1]
    beta <- res$param[2]
    Q <- informative_masking_by_rank(ts,alpha,beta)
    KL <- kl_divergence_bernoulli(P,Q)
    return(list(Q=Q,P=P,KL=KL,alpha=alpha,beta=beta,ts=ts,k=k,converged=TRUE))
}
