library(algebraic.mle)
library(stats)
library(md.tools)
library(dplyr)
library(ggplot2)

theta <- c(1.0, 1.5, 0.8, 1.25, 0.75, 1.35, 2.3)
m <- length(theta)

#set.seed(7231) # set seed for reproducibility
n <- 80
comp_times <- matrix(nrow=n,ncol=m)
for (j in 1:m)
    comp_times[,j] <- rexp(n,theta[j])
comp_times <- md_encode_matrix(comp_times,"t")

q <- 0.25
tau <- rep(-(1/sum(theta))*log(q),n)
data <- comp_times %>% md_series_lifetime_right_censoring(tau)
print(data,n=4)

############# exact componenent cause model ################
data.best <- comp_times %>%
    md_series_lifetime_right_censoring(tau) %>%
    md_bernoulli_cand_C1_C2_C3(p=0) %>%
    md_cand_sampler()

print(md_boolean_matrix_to_charsets(data.best,drop_set=TRUE),drop_latent=TRUE,n=50)

theta0 <- runif(m, min=.1, max=10)
ll.best <- md_loglike_exp_series_C1_C2_C3(data.best)
score.best <- md_score_exp_series_C1_C2_C3(data.best)
fim.best <- md_info_exp_series_C1_C2_C3(data.best)

#start <- optim(
#  par=theta0,
#  fn=ll.best,
#  method="SANN",
#  control=list(maxit=1000000,fnscale=-1,temp=100,trace=1,REPORT=1000))

start.best <- sim_anneal(
    f=ll.best,
    x0=theta0,
    options=list(
        t_init=100,
        t_end=.00001,
        alpha=.975,
        iter_per_temp=200,
        sup=function(theta) all(theta > 0),
        trace=TRUE))

loglikes.best <- apply(start.best$path,1,ll.best)
plot(loglikes.best,
    xlab="Iteration",
    ylab="Log-likelihood",
    main="Simulated Annealing: Log-likelihood (best) vs Iteration")

#mle.best.optim <- mle_optim(optim(
#    par = start.best$argmax,
#    fn = ll.best,
#    gr = score.best,
#    method = "BFGS",
#    control = list(fnscale = -1, maxit = 1000, trace = 1, REPORT=1),
#    hessian = TRUE))

mle.best.nr <- mle_newton_raphson(
    ll = ll.best,
    theta0 = start.best$argmax,
    info = fim.best,
    score = score.best,
    options = list(
        max_iter = 100000,
        eta = .05,
        trace = TRUE,
        rel_tol = 1e-9))

mle.best.nr$path
mle.best.nr$path <- mle.best.nr$path[1:44,]
loglikes.best.nr <- apply(mle.best.nr$path,1,ll.best)
plot(loglikes.best.nr,
    type="l",
    xlab="Iteration",
    ylab="Log-likelihood",
    main="Newton-Raphson: Log-likelihood (best) vs Iteration")

############## no candidate set model ################
data.worst <- comp_times %>%
    md_series_lifetime_right_censoring(tau) %>%
    md_bernoulli_cand_C1_C2_C3(p=1) %>%
    md_cand_sampler()

print(md_boolean_matrix_to_charsets(data.worst,drop_set=TRUE),drop_latent=TRUE,n=50)

ll.worst <- md_loglike_exp_series_C1_C2_C3(data.worst)
score.worst <- md_score_exp_series_C1_C2_C3(data.worst)
fim.worst <- md_info_exp_series_C1_C2_C3(data.worst)

start.worst <- sim_anneal(
    f=ll.worst,
    x0=theta0,
    options=list(
        t_init=100,
        t_end=.00001,
        alpha=.975,
        iter_per_temp=200,
        sup=function(theta) all(theta > 0),
        trace=TRUE))

loglikes.worst <- apply(start.worst$path,1,ll.worst)
plot(loglikes.worst,
    xlab="Iteration",
    ylab="Log-likelihood",
    main="Simulated Annealing: Log-likelihood (worst) vs Iteration")

mle.worst.nr <- mle_newton_raphson(
    ll = ll.worst,
    theta0 = start.worst$argmax,
    info = fim.worst,
    score = score.worst,
    options = list(
        max_iter = 100000,
        eta = 1,
        trace = TRUE,
        rel_tol = 1e-4))
    
loglikes.worst.nr <- apply(mle.worst.nr$path,1,ll.worst)
plot(loglikes.worst.nr,
    xlab="Iteration",
    ylab="Log-likelihood",
    main="Newton-Raphson: Log-likelihood (worst) vs Iteration")

############## candidate set model ################
p <- .25
data.cand <- comp_times %>%
    md_series_lifetime_right_censoring(tau) %>%
    md_bernoulli_cand_C1_C2_C3(p=p) %>%
    md_cand_sampler()

print(md_boolean_matrix_to_charsets(data.cand,drop_set=TRUE),drop_latent=TRUE,n=50)

ll.cand <- md_loglike_exp_series_C1_C2_C3(data.cand)
score.cand <- md_score_exp_series_C1_C2_C3(data.cand)
fim.cand <- md_info_exp_series_C1_C2_C3(data.cand)

start.cand <- sim_anneal(
    f=ll.cand,
    x0=theta0,
    options=list(
        t_init=100,
        t_end=.00001,
        alpha=.975,
        iter_per_temp=200,
        sup=function(theta) all(theta > 0),
        trace=TRUE))
loglikes.cand <- apply(start.cand$path,1,ll.cand)
plot(loglikes.cand,
    xlab="Iteration",
    ylab="Log-likelihood",
    main="Simulated Annealing: Log-likelihood (cand) vs Iteration")

mle.cand.nr <- mle_newton_raphson(
    ll = ll.cand,
    theta0 = start.cand$argmax,
    info = fim.cand,
    score = score.cand,
    options = list(
        max_iter = 100000,
        eta = .01,
        trace = TRUE,
        rel_tol = 1e-8))
loglikes.cand.nr <- apply(mle.cand.nr$path,1,ll.cand)
plot(loglikes.cand.nr,
    xlab="Iteration",
    ylab="Log-likelihood",
    main="Newton-Raphson: Log-likelihood (cand) vs Iteration")

### comparisons
confint(mle.best.nr)
confint(mle.worst.nr)
confint(mle.cand.nr)
sum(abs(point(mle.cand.nr)-theta))
sum(abs(point(mle.worst.nr)-theta))
sum(abs(point(mle.best.nr)-theta))
theta


