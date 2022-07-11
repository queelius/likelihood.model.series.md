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

# sim parameters
set.seed(172345)
N <- 100000
theta <- c(1,1.25,1.75,0.75,2.0,1.5,0.5)
m <- length(theta)

# we want to set delta to occur roughly 25% of the time, so we solve F(t) = .8
q.20 <- -log(1-.80)/sum(theta)

frob.norm <- function(x) sqrt(tr(t(x) %*% x))

stats2 <- tibble(n = numeric(),
                 mse = numeric(),
                 frob = numeric(),
                 rate1.lb = numeric(), rate1.ub = numeric(), rate1.in = numeric(),
                 rate2.lb = numeric(), rate2.ub = numeric(), rate2.in = numeric(),
                 rate3.lb = numeric(), rate3.ub = numeric(), rate3.in = numeric(),
                 rate4.lb = numeric(), rate4.ub = numeric(), rate4.in = numeric(),
                 rate5.lb = numeric(), rate5.ub = numeric(), rate5.in = numeric(),
                 rate6.lb = numeric(), rate6.ub = numeric(), rate6.in = numeric(),
                 rate7.lb = numeric(), rate7.ub = numeric(), rate7.in = numeric())

for (n in seq(from=100,to=N,by=100))
{
    # generate simulated masked data
    md <- tibble(t1=stats::rexp(n,rate=theta[1]),
                 t2=stats::rexp(n,rate=theta[2]),
                 t3=stats::rexp(n,rate=theta[3]),
                 t4=stats::rexp(n,rate=theta[4]),
                 t5=stats::rexp(n,rate=theta[5]),
                 t6=stats::rexp(n,rate=theta[6]),
                 t7=stats::rexp(n,rate=theta[7])) %>%
        md_series_lifetime() %>%
        md_series_lifetime_right_censoring(tau=rep(q.20,n)) %>%
        md_bernoulli_candidate_C1_C2_C3(m,p=function(n) rep(.25,n))

    l <- md_loglike_exp_series_C1_C2_C3(md)
    mle <- mle_gradient_ascent(l=l,theta0=theta)

    cat("n =",n,"mle =",t(point(mle)),"\n")

    ## bootstrap
    #mle.boot <- mle_boot_loglike(mle=mle,
    #                             loglike.gen=md_loglike_exp_series_C1_C2_C3,
    #                             data=md,
    #                             R=N,
    #                             method=mle_gradient_ascent)

    stats2 <- stats2 %>% add_row(n=n,
                                 mse=mse(mle),
                                 frob=frob.norm(vcov(mle)-ginv(md_info_exp_series_C1_C2_C3(md)(theta))),
                                 rate1.lb=confint(mle)[1,1],
                                 rate1.ub=confint(mle)[1,2],
                                 rate1.in=theta[1] >= rate1.lb && theta[1] <= rate1.ub,
                                 rate2.lb=confint(mle)[2,1],
                                 rate2.ub=confint(mle)[2,2],
                                 rate2.in=theta[2] >= rate2.lb && theta[2] <= rate2.ub,
                                 rate3.lb=confint(mle)[3,1],
                                 rate3.ub=confint(mle)[3,2],
                                 rate3.in=theta[3] >= rate3.lb && theta[3] <= rate3.ub,
                                 rate4.lb=confint(mle)[4,1],
                                 rate4.ub=confint(mle)[4,2],
                                 rate4.in=theta[4] >= rate4.lb && theta[4] <= rate4.ub,
                                 rate5.lb=confint(mle)[5,1],
                                 rate5.ub=confint(mle)[5,2],
                                 rate5.in=theta[5] >= rate5.lb && theta[5] <= rate5.ub,
                                 rate6.lb=confint(mle)[6,1],
                                 rate6.ub=confint(mle)[6,2],
                                 rate6.in=theta[6] >= rate6.lb && theta[6] <= rate6.ub,
                                 rate7.lb=confint(mle)[7,1],
                                 rate7.ub=confint(mle)[7,2],
                                 rate7.in=theta[7] >= rate7.lb && theta[7] <= rate7.ub)

    #stats2 <- stats2 %>% add_row(n=n,
    #                             asymptotic.mse=mse(mle),
    #                             boot.mse=mse(mle.boot))
}

write_csv(stats2, "data-raw/exp_series_stats_2.csv")
usethis::use_data(stats2, overwrite=T)

