rexp_series_mle <- function(n,x) {
    N <- 100000
    rate.hat <- algebraic.mle::point(x)
    sigma <- vcov(x)
    rates <- rmvnorm(N, rate.hat, sigma)
    t <- numeric(length=N)
    for (i in 1:N)
    {
        t[i] <- rexp(1, rates[i,])
    }
    t
}

N <- 1000000
t1 <- rexp_series_mle(N,rate.hat)
t1 <- t1[t1 < 1.5]
t2 <- rexp_series(N,point(rate.hat))
png("exp_series_histo")
hist(t1, main = "Histograms: ", xlab = "T ~ f(t)", col = rgb(0,0,1,0.5), ylim=c(0,5),xlim = c(0, 1.5), freq=F)
hist(t2, col = rgb(1,0,0,0.5), add = TRUE, freq=F)
dev.off()
