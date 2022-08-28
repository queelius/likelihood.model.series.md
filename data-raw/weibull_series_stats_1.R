library(dplyr)
library(tibble)
library(algebraic.mle)
library(series.system.estimation.masked.data)
library(boot)
library(usethis)
library(stats)
library(readr)
library(matlib)
library(MASS)
library(rbenchmark)

# sim parameters
set.seed(3725345)
N <- 10000
theta <- c(1.5,  10.25, # scale=lambda, shape=kappa
           1.75, 9.75,
           1.8,  10,
           1.5,  17.5)
m <- 4

 weibull_series_stats_1 <- tibble(
     n = numeric(),
     mse = numeric(),
     frob = numeric(),
     rate1.lb = numeric(), rate1.ub = numeric(), rate1.in = numeric(),
     kap1.lb  = numeric(), kap1.ub = numeric(),  kap1.in = numeric(),
     rate2.lb = numeric(), rate2.ub = numeric(), rate2.in = numeric(),
     kap2.lb  = numeric(), kap2.ub = numeric(),  kap2.in = numeric(),
     rate3.lb = numeric(), rate3.ub = numeric(), rate3.in = numeric(),
     kap3.lb  = numeric(), kap3.ub = numeric(),  kap3.in = numeric(),
     rate4.lb = numeric(), rate4.ub = numeric(), rate4.in = numeric(),
     kap4.lb  = numeric(), kap4.ub = numeric(),  kap4.in = numeric())

 for (n in seq(from=100,to=N,by=100))
{
    # generate simulated masked data
    md <- tibble(t1=stats::rweibull(n,scale=theta[1],shape=theta[2]),
                 t2=stats::rweibull(n,scale=theta[3],shape=theta[4]),
                 t3=stats::rweibull(n,scale=theta[5],shape=theta[6]),
                 t4=stats::rweibull(n,scale=theta[7],shape=theta[8])) %>%
        md_series_lifetime() %>%
        md_series_lifetime_right_censoring(tau=rep(1.0,n)) %>%
        md_bernoulli_candidate_C1_C2_C3(m,p=function(n) rep(.3,n))

    l <- md_loglike_weibull_series_C1_C2_C3(md)
    s <- md_loglike_weibull_series_score(md)
    l3 <- md_loglike_weibull_series_C1_C2_C3_3(md)

    n <- 1000
    benchmark(
        #mle <- mle_gradient_ascent(l=l,theta0=theta)
        mle <- mle_newton_raphson(l=l,theta0=theta,debug=T)
    )

    benchmark(
        mle3 <- mle_newton_raphson(l=l3,theta0=theta,debug=T)
        #mle3 <- mle_gradient_ascent(l=l3,theta0=theta)
    )


    cat("n =",n,"mle =",t(point(mle)),"\n")
    cat("n =",n,"mle2 =",t(point(mle2)),"\n")

    weibull_series_stats_1 <- weibull_series_stats_1 %>% add_row(
        n=n,
        mse=mse(mle),

        kap1.lb=confint(mle)[1,1], kap1.ub=confint(mle)[1,2],
        kap1.in=theta[1] >= kap1.lb && theta[1] <= kap1.ub,
        rate1.lb=confint(mle)[2,1], rate1.ub=confint(mle)[2,2],
        rate1.in=theta[2] >= rate1.lb && theta[2] <= rate1.ub,

        kap2.lb=confint(mle)[3,1], kap2.ub=confint(mle)[3,2],
        kap2.in=theta[3] >= kap1.lb && theta[3] <= kap2.ub,
        rate2.lb=confint(mle)[4,1], rate2.ub=confint(mle)[4,2],
        rate2.in=theta[4] >= rate2.lb && theta[4] <= rate2.ub,

        kap3.lb=confint(mle)[5,1], kap3.ub=confint(mle)[5,2],
        kap3.in=theta[5] >= kap3.lb && theta[5] <= kap3.ub,
        rate3.lb=confint(mle)[6,1], rate3.ub=confint(mle)[6,2],
        rate3.in=theta[6] >= rate3.lb && theta[6] <= rate3.ub,

        kap4.lb=confint(mle)[7,1], kap4.ub=confint(mle)[7,2],
        kap4.in=theta[7] >= kap4.lb && theta[7] <= kap4.ub,
        rate4.lb=confint(mle)[8,1], rate1.ub=confint(mle)[8,2],
        rate4.in=theta[8] >= rate4.lb && theta[8] <= rate4.ub)
}

write_csv(weibull_series_stats_2, "data-raw/weibull_series_stats_2.csv")
usethis::use_data(weibull_series_stats_2, overwrite=T)

