library(dplyr)
library(tibble)

##### simulation parameters ######
wei.true <- c(3,4,
              5,6,
              7,8) # set the true parameter values
n <- 38               # set the sample size
N <- 1000             # set the number of simulations
tau <- .22             # right-censoring time
p <- 0.43              # probability a non-failed component is in candidate set

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

# Create an empty vector to store the estimated parameter values
theta.hats <- matrix(nrow=N,ncol=3)
theta.mles <- list(length=N)
# Simulate n_sim samples from the normal distribution with the true parameter values
for (i in 1:N)
{
    md <- data.frame(
        t1=rweibull(n,wei.true[1],wei.true[2]),
        t2=rweibull(n,wei.true[3],wei.true[4]),
        t3=rweibull(n,wei.true[5],wei.true[6])) %>%
        md_series_lifetime_right_censoring(tau) %>%
        md_bernoulli_cand_C1_C2_C3(p) %>%
        md_cand_sampler()

    # find the MLE of rate parameters
    # we help MLE method by initially guessing the true parameter
    theta.hat <- md_mle_weibull_series_C1_C2_C3(md,wei.true)
    theta.hats[i,] <- point(theta.hat)
    theta.mles[[i]] <- theta.hat
    print(theta.hats[i,])
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
png("wei_series_samp_dist")
scatterplot3d(theta.hats[,1],
              theta.hats[,2],
              theta.hats[,3],
              main = expression(
                  paste("Sampling distribution of ", hat(theta)))
dev.off()

# Calculate the bias of the MLE estimator
bias <- colMeans(theta.hats - wei.true)

# MSE
mse <- colMeans((theta.hats - wei.true)^2)

cat("Bias of the MLE estimator:", bias)
cat("MSE component wise:", mse)
cat("total MSE:", mean(mse))
cat("MSE from asymptotic theory:", mse(theta.hat))

# compute CI coverage prob
cis <- list(length=N)
for (i in 1:N)
{
    cis <- confint(theta.mles[[i]])
}
