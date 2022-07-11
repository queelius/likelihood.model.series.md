library(dplyr)
library(tibble)
library(algebraic.mle)
#library(series.system.estimation.masked.data)
library(boot)
library(usethis)
library(stats)
library(readr)

# sim parameters
set.seed(1234)
N <- 5000
theta <- c(1,1.25,1.75)
m <- length(theta)

# we want to set delta to occur roughly 25% of the time, so we solve F(t) = .8
q.25 <- -log(1-.75)/sum(theta)

exp_series_stats_1 <- tibble(
    n = numeric(),
    asymptotic.mse = numeric(), boot.mse = numeric(),
    asymptotic.rate1.se = numeric(),  boot.rate1.se = numeric(),
    asymptotic.rate2.se = numeric(),  boot.rate2.se = numeric(),
    asymptotic.rate3.se = numeric(),  boot.rate3.se = numeric(),
    rate1 = numeric(), rate1.bias = numeric(),
    rate2 = numeric(), rate2.bias = numeric(),
    rate3 = numeric(), rate3.bias = numeric())

for (n in seq(from=100,to=N,by=100))
{
    # generate simulated masked data
    md <- tibble(t1=stats::rexp(n,rate=theta[1]),
                 t2=stats::rexp(n,rate=theta[2]),
                 t3=stats::rexp(n,rate=theta[3])) %>%
        md_series_lifetime() %>%
        md_series_lifetime_right_censoring(tau=rep(q.25,n)) %>%
        md_bernoulli_candidate_C1_C2_C3(m,p=function(n) rep(.3,n))

    l <- md_loglike_exp_series_C1_C2_C3(md)
    scr <- md_score_exp_series_C1_C2_C3(md)
    mle <- mle_gradient_ascent(l=l,theta0=theta,score=scr)

    cat("n =",n,"mle =",t(point(mle)),"\n")

    ## bootstrap
    mle.boot <- mle_boot_loglike(mle=mle,
                                 loglike.gen=md_loglike_exp_series_C1_C2_C3,
                                 data=md,
                                 R=N)

    exp_series_stats_1 <- exp_series_stats_1 %>% add_row(n=n,
                               asymptotic.mse=mse(mle), boot.mse=mse(mle.boot),
                               asymptotic.rate1.se=se(mle)[1], boot.rate1.se=se(mle.boot)[1],
                               asymptotic.rate2.se=se(mle)[2], boot.rate2.se=se(mle.boot)[2],
                               asymptotic.rate3.se=se(mle)[3], boot.rate3.se=se(mle.boot)[3],
                               rate1=point(mle)[1], rate1.bias=bias(mle.boot)[1],
                               rate2=point(mle)[2], rate2.bias=bias(mle.boot)[2],
                               rate3=point(mle)[3], rate3.bias=bias(mle.boot)[3])

    print(exp_series_stats_1[nrow(exp_series_stats_1),])
}

write_csv(exp_series_stats_1, "data-raw/exp_series_stats_1.csv")
usethis::use_data(exp_series_stats_1, overwrite=T)

