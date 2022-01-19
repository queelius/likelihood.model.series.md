# series system with m=3 nodes, each lomax distributed
# system parameter value: lambda=(3,4,5), kappa=(2,3,4)
# sample size: n=10000
# candidate model m0
# each candidate set has w=2 candidates
# starting seed is 5123, so it will always produce the same output

library(readr)
library(masked.data)
library(usethis)

lambda <- c(3,4,5)
kappa <- c(2,3,4)
n <- 10000
w <- rep(2,n)
set.seed(5123)

lomax_series_data_1 <- masked.data::md_lomax_series(
    n=n,
    lambda=lambda,
    kappa=kappa,
    w=w,
    candidate_model=md_candidate_m0,
    metadata=T)

masked.data::md_write_csv(lomax_series_data_1, "data-raw/lomax_series_data_1.csv")
usethis::use_data(lomax_series_data_1, overwrite = TRUE)
