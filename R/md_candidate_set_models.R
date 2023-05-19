#' md_bernoulli_cand_C1_kld
#'
#' For each row (observation) in `md`, a probability `Q = (q1, q2, ..., qm)` is
#' constructed such that `qj` is the probability that the j-th component is in
#' the candidate set, `qk = 1`, where `k` is failed component.
#'
#' `Q` is an *informed* candidate model that uses `informative_masking_by_rank`
#' to assign higher probabilities to components that failed earlier (which is
#' something we typically only know in, say, a simulation study).
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
#' @param md component failure times for the series system
#' @param p numeric, defines `P = (p, ..., p, 1, p, ..., p)`.
#' @param d numeric, the KL divergence from P = (p, p, ..., p, 1, p, ..., p)
#'          to try to obtain
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
#' @importFrom md.tools md_decode_matrix
#' @export
md_bernoulli_cand_C1_kld <- function(
        md, p, d, eps=1e-4, max_iter=100000L, lr=1, lambda=1, alpha0=5, beta0=.5, debug=F)
{
    t <- md_decode_matrix(md,"t")
    m <- ncol(t)
    n <- nrow(t)
    stopifnot(n > 0, m > 0)

    for (i in 1:nrow(md))
    {
        if (debug)
            cat("trial ", i, "\n")

        row <- md[i,]
        tryCatch(expr = {
            res <- generate_relaxed_cand_C1(d=d,
                                            ts=t[i,],
                                            p=p,
                                            debug=debug,
                                            eps=eps,
                                            lr=lr,
                                            max_iter=max_iter,
                                            lambda=lambda,
                                            alpha0=alpha0,
                                            beta0=beta0)
            if (res$converged && abs(sum(res$Q) - sum(res$P)) < .05)
            {
                kl.d <- kl_divergence_bernoulli(res$P,res$Q)
                md[i, "d"] <- kl.d
                md[i, "p"] <- p
                md[i, paste0("q",1:m)] <- t(res$Q)

                if (abs(kl.d - d) > 1e-2)
                {
                    cat("Violation of `KL - d` soft constraint: ",
                        kl.d - d, "\n")
                }
            }
            else
            {
                cat("Failed to converge: stopped at (alpha,beta) = (",
                    res$alpha, ", ", res$beta, ")\n")
            }
        },
        error = function(e) {
            cat("Caught an error: ", e$message, "\n")
        })
    }
    md
}

#' md_bernoulli_cand_C1_C2_C3
#'
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
#' @param p a parameter vector of probabilities
#' @param compvar column name of the component lifetime variables, defaults to
#'                `t`, e.g., `t1, t2, ..., tm`.
#' @param qvar column prefix for component probabilities, defaults to `q`,
#'             e.g., `q1, q2, ..., qm`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @export
md_bernoulli_cand_C1_C2_C3 <- function(md,p,compvar="t",qvar="q")
{
    Tm <- md_decode_matrix(md, compvar)
    stopifnot(!is.null(Tm))
    m <- ncol(Tm)
    n <- nrow(Tm)
    stopifnot(n > 0, m > 0)

    Q <- matrix(rep(p,n*m), nrow=n)
    #for (i in 1:n)
    #    Q[i,which.min(Tm[i,])] <- 1
    Q[cbind(1:n, apply(Tm, 1, which.min))] <- 1
    if ("delta" %in% colnames(md))
    {
        Q[which(md$delta),] <- 0
        #for (i in 1:n)
        #{
        #    if (md$delta[i])
        #        Q[i,] <- rep(0,m)
        #}
    }

    md[,paste0(qvar,1:m)] <- NULL
    md %>% bind_cols(md_encode_matrix(Q,qvar)) %>% md_mark_latent(paste0(qvar,1:m))
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

#' md_block_candidate_m3
#'
#' Block candidate model. It produces an MLE that is non-unique, i.e.,
#' if estimating exponential series system, the MLE
#' for λ₁, λ₂, and λ₃ satisfies
#'     `λ̂₁ + λ̂₂ = λ₁ + λ₂`
#' and
#'     `λ̂₃ = λ₃`.
#'
#' This illustrates one of the conditions for the MLE to be unique:
#'
#'     (1) multiple components must not always show together. if they
#'     do, then the best we can learn is that they satisfy some hyper-plane
#'     constraint.
#'
#' @param md masked data
#' @param compvar Prefix-code for the component lifetimes, defaults to `t`,
#'                e.g., `t1, t2, ..., tm`.
#' @param qvar column prefix for component probabilities, defaults to `q`,
#'             e.g., `q1, q2, ..., qm`.
#' @importFrom md.tools md_decode_matrix md_encode_matrix
#' @importFrom dplyr %>% bind_cols
#' @export
md_block_candidate_m3 <- function(md,qvar="q",compvar="t")
{
    Tm <- md_decode_matrix(md,"t")
    stopifnot(!is.null(Tm))
    m <- ncol(Tm)
    n <- nrow(Tm)
    stopifnot(m == 3, n > 0)

    block <- function(ts)
    {
        k <- which.min(ts)
        if (k == 1)
            return(c(1,1,0))
        if (k == 2)
            return(c(1,1,0))
        if (k == 3)
        {
            if (runif(1) < 0.1)
                return(c(1,1,1))
            else
                return(c(0,0,1))
        }
    }


    Q <- matrix(rep(NA,m*n), nrow=n)
    for (i in 1:n)
        Q[i,] <- block(Tm[i,])

    md %>% bind_cols(md_encode_matrix(Q,qvar))
}

