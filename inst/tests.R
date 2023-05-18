library(dplyr)
library(tibble)
library(algebraic.mle)
library(stats)
library(microbenchmark)

# sim parameters
set.seed(5213334)
theta <- c(1,1.25,1.75)
m <- length(theta)
n <- 100000
# we want to set delta to occur roughly alpha=25% of the time,
# so we solve for t in the expression F(t) = alpha
alpha <- .25
q.alpha <- -log(1-alpha)/sum(theta)

# with probability p=.5, a non-failed component
# appears in the candidate set (bernoulli model for
# simulation purposes, which satisfies regularity
# conditions).
p <- rep(.5,n)

# generate simulated masked data according to
# simulation parameters
md <- tibble(t1=stats::rexp(n,rate=theta[1]),
             t2=stats::rexp(n,rate=theta[2]),
             t3=stats::rexp(n,rate=theta[3])) %>%
    md_series_lifetime() %>%
    md_series_lifetime_right_censoring(tau=rep(q.alpha,n)) %>%
    md_bernoulli_candidate_C1_C2_C3(m,p)

loglik <- md_loglike_exp_series_C1_C2_C3(md)
benchmark_results <- microbenchmark(
    ,
    ,
    times = 10)

# Print the results
print(benchmark_results)

(mle.fast <- md_mle_exp_series_C1_C2_C3(md=md,theta0=theta,eps=1e-30))
(mle.slow <- md_mle_exp_series_C1_C2_C3(md=md,theta0=theta,eps=1e-5,method="slow"))
(mle1 <- mle_constrained_solver(nll = function(theta) -loglik(theta),
                                theta0 = theta0,
                                eps=1e-10))


(l1 <- loglik(point(mle2)))
(l2 <- loglik(mle1$par))

md_mle_exp_series_C1_C2_C3
