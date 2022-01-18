# series system with m=5 nodes, each exponentially distributed
# system parameter value: rate.star=(3,5,4,6,7)
# sample size: n=1000
# candidate model m0
# candidate set sizes are randomly drawn from {2,3,4}
# starting seed is 12345, so it will always produce the same output

library(readr)
library(masked.data)
library(usethis)

set.seed(12345)
rate.star <- c(3,5,4,6,7)
n <- 100000
w <- sample(2:4,n,replace=T)

data2 <- masked.data::md_exp_series(
    n=n,
    theta=rate.star,
    w=w,
    candidate_model=md_candidate_m0,
    metadata=T)

masked.data::md_write_csv(data2, "data-raw/data2.csv")
usethis::use_data(data2, overwrite = TRUE)
