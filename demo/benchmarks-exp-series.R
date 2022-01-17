library(rbenchmark)
library(masked.data)
library(ggplot2)

n <- 100
w <- rep(2,n)
md <- masked.data::md_exp_series(n=n,theta=rate.star,w=w)
ll.m0 <- masked.data::md_kloglike_exp_series_m0(md)
ll.m0.ref <- masked.data::md_kloglike_exp_series_m0_ref(md)

res <- microbenchmark(
    "1" =
    {
        optim(runif(3,.01,10),function(rate) { -ll.m0.ref(rate) })
    },
    "2" =
    {
        optim(runif(3,.01,10),function(rate) { -ll.m0(rate)})
    })

print(res)
autoplot(res)
