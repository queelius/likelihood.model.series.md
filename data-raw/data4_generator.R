# series system with m=3 nodes, each exponentially distributed
# system parameter value: rate.star=(3,3,3)
# sample size: n=10000
# candidate model m0
# candidate set sizes are randomly drawn from {2}
# starting seed is 345, so it will always produce the same output

library(readr)
library(masked.data)
library(usethis)

set.seed(345)
rate.star <- c(1,1,1)
n <- 10000
w <- rep(2,n)

data4 <- masked.data::md_exp_series(
    n=n,
    theta=rate.star,
    w=w,
    candidate_model=md_candidate_m0,
    metadata=T)

masked.data::md_write_csv(data2, "data-raw/data4.csv")
usethis::use_data(data4, overwrite = TRUE)
