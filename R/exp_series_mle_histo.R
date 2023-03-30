#' suppose (T,rate.hat) ~ f(t,r). then, we can find the marginal of T
#' by integrating over rate.hat.
#' we just do Monte Carlo simulation, sampling from joint (T,rate.hat)
#' and taking only the observed T values
rexp_series_mle <- function(n,x) {
    N <- 100000
    rate.hat <- algebraic.mle::point(x)
    sigma <- vcov(x)
    rates <- rmvnorm(N, rate.hat, sigma)
    t <- numeric(length=N)
    for (i in 1:N)
        t[i] <- rexp(1, rates[i,])
    t
}

# sample size
N <- 1000000

# sample from T
t1 <- rexp_series_mle(N,rate.hat)
t1 <- t1[t1 < 1.5]

# sample from T | rate.hat ~ EXP(rate.hat)
t2 <- rexp_series(N,point(rate.hat))

# plot histo grams
# and save to file `exp_series_histo.png`
png("exp_series_histo")
hist(t1,
     main = expression(
         paste("Histograms: T ~ f(t) vs T|",
                   hat(lambda), " ~ f(t|",
                   hat(lambda), ")")),
     xlab = "T", col = rgb(0,0,1,0.5),
     ylim=c(0,5),xlim = c(0, 1.5), freq=F)
hist(t2, col = rgb(1,0,0,0.5), add = TRUE, freq=F)
dev.off()
