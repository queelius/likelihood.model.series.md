library(usethis)
library(stats)
library(dplyr)
library(readr)

# sim parameters
set.seed(123)
n <- 200
m <- 3
rates <- c(3,4,5)

# we want to set delta to occur roughly 20% of the time, so we solve F(t) = .8
q.20 <- -log(1-.80)/sum(rates)
tau <- rep(q.20,n)
p <- runif(n)

# generate simulated masked data
exp_series_md_1 <- tibble(t1=stats::rexp(n,rate=rates[1]),
                          t2=stats::rexp(n,rate=rates[2]),
                          t3=stats::rexp(n,rate=rates[3])) %>%
    md_series_lifetime() %>%
    #md_series_lifetime_right_censoring(tau) %>%
    md_bernoulli_candidate_C1_C2_C3(m,p)

exp_series_md_1$tau <- NULL
write_csv(exp_series_md_1, "data-raw/exp_series_md_1.csv")
usethis::use_data(exp_series_md_1, overwrite=T)

