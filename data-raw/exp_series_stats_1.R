library(dplyr)
library(tibble)
library(algebraic.mle)
#library(series.system.estimation.masked.data)
library(boot)
library(usethis)
library(stats)
library(readr)

# sim parameters
set.seed(12334)
N <- 50000
theta <- c(1,1.25,1.75)
m <- length(theta)
step_size = 100
# we want to set delta to occur roughly 25% of the time, so we solve F(t) = .8
q.25 <- -log(1-.75)/sum(theta)

exp_series_stats_1 <- tibble(
    N = numeric()/step_size,
    asymptotic.mse = numeric(), boot.mse = numeric(),
    asymptotic.rate1.se = numeric(),  boot.rate1.se = numeric(),
    asymptotic.rate2.se = numeric(),  boot.rate2.se = numeric(),
    asymptotic.rate3.se = numeric(),  boot.rate3.se = numeric(),
    rate1 = numeric(), rate1.bias = numeric(),
    rate2 = numeric(), rate2.bias = numeric(),
    rate3 = numeric(), rate3.bias = numeric())

for (n in seq(from=step_size,to=N,by=step_size))
{
    p <- rep(.333,n)
    # generate simulated masked data
    md <- tibble(t1=stats::rexp(n,rate=theta[1]),
                 t2=stats::rexp(n,rate=theta[2]),
                 t3=stats::rexp(n,rate=theta[3])) %>%
        md_series_lifetime() %>%
        md_series_lifetime_right_censoring(tau=rep(q.25,n)) %>%
        md_bernoulli_candidate_C1_C2_C3(m,p)

    mle <- md_mle_exp_series_C1_C2_C3(md=md,theta0=theta)
    cat("n =",n,"mle =",t(point(mle)),"\n")

    solver <- function(md,index,...) {
        theta.hat <- md_mle_exp_series_C1_C2_C3(md[index,],theta,...)
        count <<- count + 1
        return(theta.hat)
    }

    ## bootstrap
    mle.boot <- boot::boot(
        statistic=solver,
        data=md,
        R=999)

    #
    # exp_series_stats_1 <- exp_series_stats_1 %>% add_row(n=n,
    #                            asymptotic.mse=mse(mle), boot.mse=mse(mle.boot),
    #                            asymptotic.rate1.se=se(mle)[1], boot.rate1.se=se(mle.boot)[1],
    #                            asymptotic.rate2.se=se(mle)[2], boot.rate2.se=se(mle.boot)[2],
    #                            asymptotic.rate3.se=se(mle)[3], boot.rate3.se=se(mle.boot)[3],
    #                            rate1=point(mle)[1], rate1.bias=bias(mle.boot)[1],
    #                            rate2=point(mle)[2], rate2.bias=bias(mle.boot)[2],
    #                            rate3=point(mle)[3], rate3.bias=bias(mle.boot)[3])
    #
    # print(exp_series_stats_1[nrow(exp_series_stats_1),])
}

write_csv(exp_series_stats_1, "data-raw/exp_series_stats_1.csv")
usethis::use_data(exp_series_stats_1, overwrite=T)

n <- 1000
p <- rep(.1,n)
# generate simulated masked data
md <- tibble(t1=stats::rexp(n,rate=theta[1]),
             t2=stats::rexp(n,rate=theta[2]),
             t3=stats::rexp(n,rate=theta[3])) %>%
    md_series_lifetime() %>%
    md_series_lifetime_right_censoring(tau=rep(q.25,n)) %>%
    md_bernoulli_candidate_C1_C2_C3(m,p)
md_mle_exp_series_C1_C2_C3(md,theta,method="slow")
md_mle_exp_series_C1_C2_C3(md,theta)

solver <- function(theta,ind,...)
{
    point(md_mle_exp_series_C1_C2_C3(md[ind,],theta,method="slow"))
}

mle.boot <- boot::boot(statistic=solver,data=md,R=999)

