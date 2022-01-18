# series system with m=3 nodes, each exponentially distributed
# system parameter value: rate.star=(3,4,5)
# sample size: n=1000
# candidate model m0
# each candidate set has w=2 candidates
# starting seed is 123, so it will always produce the same output

library(readr)
library(masked.data)
library(usethis)

rate.star <- c(3,4,5)
n <- 1000
w <- rep(2,n)
set.seed(123)

data1 <- masked.data::md_exp_series(
    n=n,
    theta=rate.star,
    w=w,
    candidate_model=md_candidate_m0,
    metadata=T)

masked.data::md_write_csv(data1, "data-raw/data1.csv")
usethis::use_data(data1, overwrite = TRUE)
