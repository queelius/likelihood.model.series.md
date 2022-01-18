library(microbenchmark)
library(masked.data)
library(ggplot2)

ll.m0 <- masked.data::md_kloglike_exp_series_m0(masked.data::data3)

score.m0.slow <- function(rate) { numDeriv::grad(ll.m0,rate) }
score.m0 <- md_score_exp_series_m0(data3)

info.m0.slow <- function(x) { numDeriv::hessian(function(rate) { -ll.m0(rate) },x) }
#info.m0 <- function(x) { numDeriv::jacobian(function(rate) { -score.m0(rate) },x) }
info.m0 <- md_info_exp_series_m0(masked.data::data3)

m <- md_num_nodes(masked.data::data3)


res <- microbenchmark(
#    fisher_scoring_faster =
#    {
        md_fisher_scoring(rep(1,m),info.m0.slow,score.m0,1e-3)
#    },
    fisher_scoring_fastest =
    {
        md_fisher_scoring(rep(1,m),info.m0,score.m0,1e-3)
    }
    #optim =
    #{
    #    optim(rep(1,m),function(rate) { -ll.m0(rate)})
    #},
#    fisher_scoring =
#    {
#        md_fisher_scoring(c(1,1,1,1,1),info.m0.slow,score.m0.slow,1e-3)
#    }
)

print(res)
boxplot(res)
autoplot(res)
