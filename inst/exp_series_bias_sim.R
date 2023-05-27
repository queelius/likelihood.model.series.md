library(dplyr)
library(tibble)

##### simulation parameters ######
rate.true <- c(3,4,5) # set the true parameter values
n <- 50               # set the sample size
N <- 1000             # set the number of simulations
tau <- .2             # right-censoring time
p <- 0.5              # probability a non-failed component is in candidate set

# a candidate set model that obeys the conditions 1, 2, an 3
# namely:
#     1) given ...
md_candidate_set_C1_C2_C3 <- function(md,m,p)
{
    stopifnot(!is.null(md$k))
    n <- nrow(md)
    stopifnot(n > 0)

    x <- matrix(NA,nrow=n,ncol=m)
    u <- matrix(runif(m*n),nrow=n)

    for (i in 1:n)
    {
        for (j in 1:m)
            x[i,j] <- ifelse(md$k[i]==j, T, u[i,j] < p)
    }

    x <- tibble::as_tibble(x)
    colnames(x) <- paste0("x",1:m)
    md %>% dplyr::bind_cols(x)
}

library(tidyverse)
make_exp_series_sample <- function(n,rate,tau,p)
{
    m <- length(rate)
    md <- data.frame(t1 = rexp(n,rate[1]))
    for (j in 2:m)
        md[[paste0("t",j)]] <- rexp(n,rate[j])

    md %>% md_series_lifetime() %>%
        md_series_lifetime_right_censoring(tau) %>%
        md_candidate_set_C1_C2_C3(m,p)
}
samp <- make_exp_series_sample(n,rate.true,tau,p)
md_mle_exp_series_C1_C2_C3(samp,rate.)
# Create an empty vector to store the estimated parameter values
rate.hats <- matrix(nrow=N,ncol=3)
rate.mles <- list(length=N)
# Simulate n_sim samples from the normal distribution with the true parameter values
for (i in 1:N)
{
    md <- data.frame(
        t1=rexp(n,rate.true[1]),
        t2=rexp(n,rate.true[1]),
        t3=rexp(n,rate.true[1])) %>%
        md_series_lifetime() %>%
        md_series_lifetime_right_censoring(tau) %>%
        md_candidate_set_C1_C2_C3(3,p)

    # find the MLE of rate parameters
    # we help MLE method by initially guessing the true parameter
    rate.hat <- md_mle_exp_series_C1_C2_C3(md,rate_true)
    rate.hats[i,] <- point(rate.hat)
    rate.mles[[i]] <- rate.hat
    print(rate.hats[i,])
}

# if MLE converged to a local max, this is true:
# is.positive.definite(sigma)


# do coverage probability
# bias

# mean squared error
# E((rate.hat - rate.true)'(rate.hat - rate.true))
# tr(V) - bias'bias

# variance-covariance
# E((rate.hat - rate.true)'(rate.hat - rate.true))


# plot histo grams
# and save to file `exp_series_histo.png`
png("exp_series_samp_dist")
scatterplot3d(rate.hats[,1],
              rate.hats[,2],
              rate.hats[,3],
              main = expression(
                  paste("Sampling distribution of ", hat(lambda))))
dev.off()

# Calculate the bias of the MLE estimator
bias <- colMeans(rate.hats - rate.true)

# MSE
mse <- colMeans((rate.hats - rate.true)^2)

cat("Bias of the MLE estimator:", bias)
cat("MSE component wise:", mse)
cat("total MSE:", mean(mse))
cat("MSE from asymptotic theory:", mse(rate.hat))

# compute CI coverage prob
for (i in 1:N)
{
    cis <- confint(rate.mles[[i]])
}



























