# series system with m=5 nodes, each exponentially distributed
# system parameter value: rate=(2,4,8)
# sample size: n=250
# candidate model m0
# candidate set sizes are randomly drawn from {1,2,3}
# starting seed is 1233445, so it will always produce the same output

library(readr)
library(masked.data)
library(usethis)

set.seed(1233445)
rate <- c(2,4,8)
n <- 250
w <- sample(1:3,n,replace=T)

exp_series_data_5 <- masked.data::md_exp_series(
    n=n,
    theta=rate,
    w=w,
    candidate_model=md_candidate_m0,
    metadata=T)

masked.data::md_write_csv(exp_series_data_5, "data-raw/exp_series_data_5.csv")
usethis::use_data(exp_series_data_5, overwrite = TRUE)
