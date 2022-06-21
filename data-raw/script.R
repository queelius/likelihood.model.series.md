library(dplyr)
library(tibble)
library(algebraic.mle)
library(series.system.estimation.masked.data)
library(boot)

n <- 5000
theta <- c(1,1.25,1.75)
m <- length(theta)

md <- tibble(t1=rexp(n,theta[1]),
             t2=rexp(n,theta[2]),
             t3=rexp(n,theta[3])) %>%
    md_series_lifetime() %>%
    md_bernoulli_candidate_C1_C2_C3(m, function(n) rep(.25,n))
print(md)
l <- md_loglike_exp_series_C1_C2_C3(md)
scr <- md_score_exp_series_C1_C2_C3(md)
mle <- mle_gradient_ascent(l=l,theta0=theta,score=scr)
theta.hat <- point(mle)
summary(mle)
(V.hat <- vcov(mle))
(mse <- sum(diag(V.hat)))

## bootstrap
theta.boot <- mle_boot_loglike(mle=mle,
                               loglike.gen=md_loglike_exp_series_C1_C2_C3,
                               data=data)

mle.boot <- boot(md,
                 statistic=function(x,idx) point(mle_newton_raphson(
                     l=md_loglike_exp_series_C1_C2_C3(x[idx,]),
                     theta0=theta)),
                 R=n)
