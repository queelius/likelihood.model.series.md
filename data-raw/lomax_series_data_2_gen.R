# series system with m=3 nodes, each lomax distributed
# system parameter value: lambda=(1,1.5,0.75), kappa=(2,1.5,2.5)
# sample size: n=2000
# candidate model m0
# each candidate set has w=2 candidates
# starting seed is 14523, so it will always produce the same output

library(readr)
library(masked.data)
library(usethis)

lambda <- c(1,1.5,.75)
kappa <- c(2,1.5,2.5)
n <- 2000
w <- rep(2,n)
set.seed(14523)

lomax_series_data_2 <- masked.data::md_lomax_series(
    n=n,
    lambda=lambda,
    kappa=kappa,
    w=w,
    candidate_model=md_candidate_m0,
    metadata=T)

masked.data::md_write_csv(lomax_series_data_2, "data-raw/lomax_series_data_2.csv")
usethis::use_data(lomax_series_data_2, overwrite = TRUE)
