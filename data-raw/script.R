library(QuantPsyc)
library(MVN)
library(dplyr)
library(tibble)
library(algebraic.mle)
library(series.system.estimation.masked.data)
library(boot)

n <- 5000
theta <- c(1,1.25,1.75)
m <- length(theta)

md <- tibble(t1=rexp(n,theta[1]),
             t2=rexp(n,theta[2]),
             t3=rexp(n,theta[3])) %>%
    md_series_lifetime() %>%
    md_bernoulli_candidate_C1_C2_C3(m, function(n) rep(.25,n))
print(md)
l <- md_loglike_exp_series_C1_C2_C3(md)
scr <- md_score_exp_series_C1_C2_C3(md)
mle <- mle_gradient_ascent(l=l,theta0=theta,score=scr)

## bootstrap
mle.boot <- mle_boot_loglike(mle=mle,
                             loglike.gen=md_loglike_exp_series_C1_C2_C3,
                             data=md,
                             R=n)

#mle.boot <- boot(md,
#                 statistic=function(x,idx) point(mle_newton_raphson(
#                     l=md_loglike_exp_series_C1_C2_C3(x[idx,]),
#                     theta0=theta)),
#                 R=n)

vcov(mle)
vcov(mle.boot)

mse(mle) # sum(diag(vcov(mle))))
mse(mle.boot)

se(mle)
se(mle.boot)

bias(mle)
bias(mle.boot)

point(mle)
point(mle.boot)

summary(mle)
summary(mle.boot)


sampler(mle.boot)(10)
sampler(mle)(10)


############################
theta <- c(1,1.25,1.75)
m <- length(theta)

set.seed(47937123)
N <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
mles <- matrix(nrow=length(N),ncol=m)
biases <- matrix(nrow=length(N),ncol=m)
mses <- numeric(length(N))
mses.asym <- numeric(length(N))
sds <- matrix(nrow=length(N),ncol=m)
sds.asym <- matrix(nrow=length(N),ncol=m)
covs <- list()
covs.asym <- list()
# mardia mult.norm(theta.boot$t)

for (i in 1:length(N))
{
    md <- tibble(
        t1=rexp(N[i],theta[1]),
        t2=rexp(N[i],theta[2]),
        t3=rexp(N[i],theta[3])) %>%
        md_series_lifetime() %>%
        md_bernoulli_candidate_C1_C2_C3(m=m,p=function(n) rep(.25,n))

    loglik <- md_loglike_exp_series_C1_C2_C3(md)
    mle <- mle_gradient_ascent(l=loglik,theta0=theta)
    theta.boot <- mle_boot_loglike(mle=mle,
                                   loglike.gen=md_loglike_exp_series_C1_C2_C3,
                                   data=md,
                                   R=N[i])

    mles[i,] <- as.numeric(point(mle))
    biases[i,] <- bias(x=theta.boot,par=theta)
    mses[i] <- mse(x=theta.boot,par=theta)
    sds[i,] <- sqrt(diag(vcov(theta.boot)))
    sds.asym[i,] <- sqrt(diag(vcov(mle)))
    mses.asym[i] <- mse(mle)
    covs[[i]] <- vcov(theta.boot)
    covs.asym[[i]] <- vcov(mle)



    cat("n:",         N[i],         "\n")
    cat("mle:",       mles[i,],     "\n")
    cat("sd:",        sds[i,],      "\n")
    cat("bias:",      biases[i,],   "\n")
    cat("mse:",       mses[i],      "\n")
    cat("----------------------------\n")
    cat("mse.asym:", mses.asym[i],  "\n")
    cat("sds.asym:",  sds.asym[i,], "\n")
    cat("============================\n")
}

res <- tibble(n=as.integer(N),mses=mses,bias=biases)
#readr::write_csv2(res)
#sink('boot-analysis-output.txt')





#ootstrap bca confidence intervals
#
#   2.5%   97.5 %
#   1 0.9217328 1.049846
#   1.2042913 1.348674
#   1.6523594 1.815097




mle_exp_series_solver <- function(md,theta0)
{
    l <- md_loglike_exp_series_C1_C2_C3(md)
    #s <- md_score_exp_series_C1_C2_C3(md)
    #J <- md_info_exp_series_C1_C2_C3(md)
    mle_gradient_ascent(l=l,theta0=theta0)
}

