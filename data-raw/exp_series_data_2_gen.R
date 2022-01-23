# series system with m=5 nodes, each exponentially distributed
# system parameter value: rate=(3,5,4,6,7)
# sample size: n=100000
# candidate model m0
# candidate set sizes are randomly drawn from {2,3,4}
# starting seed is 12345, so it will always produce the same output

library(readr)
library(masked.data)
library(usethis)

set.seed(12345)
rate <- c(3,5,4,6,7)
n <- 100000
w <- sample(2:4,n,replace=T)

exp_series_data_2 <- masked.data::md_exp_series(
    n=n,
    theta=rate,
    w=w,
    candidate_model=md_candidate_m0,
    metadata=T)

masked.data::md_write_csv(exp_series_data_2, "data-raw/exp_series_data_2.csv")
usethis::use_data(exp_series_data_2, overwrite = TRUE)
