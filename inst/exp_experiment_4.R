library(utils)
library(ggplot2)
library(tidyverse)
library(algebraic.mle)
library(stats)

generate_data <- function(n, theta, tau, p) {
    m <- length(theta)
    comp_times <- matrix(nrow=n,ncol=m)
    for (j in 1:m)
        comp_times[,j] <- rexp(n,theta[j])
    comp_times <- md_encode_matrix(comp_times,"t")

    comp_times %>%
        md_series_lifetime_right_censoring(tau) %>%
        md_bernoulli_cand_C1_C2_C3(p) %>%
        md_cand_sampler()
}

custom_solver <- function(data, theta, extra_info = NULL, annealing = TRUE) {
    ll <- md_loglike_exp_series_C1_C2_C3(data)
    ll.grad <- md_score_exp_series_C1_C2_C3(data)
    fish <- md_fim_exp_series_C1_C2_C3(data)
    ll.hess <- function(x) -fish(x)
    theta.hat <- NULL
    start <- NULL
    m <- length(theta)

    tryCatch({
        start <- list(par = theta)
        if (annealing) {
            start <- sim_anneal(par = theta, fn = ll, control =
                list(fnscale = -1, maxit = 100000L, trace = FALSE,
                    t_init = 20, 1e-3, alpha = 0.95,
                    it_per_temp = 50L,
                    neigh = function(par, temp, ...) {
                        tt <- min(temp, 1)
                        par + rnorm(m, 0, tt)
                    },
                    sup = function(x) all(x > 0)))
        }
        res_optim <- optim(
            par = start$par,
            fn = ll,
            gr = ll.grad,
            method = "L-BFGS-B",
            lower = 1e-30,
            hessian = FALSE,
            control = list(fnscale = -1, maxit = 1000L))
        theta.hat <- mle_numerical(res_optim,
            options = list(hessian = ll.hess(res_optim$par)))
    }, error = function(e) {
        cat("Sample size", nrow(data), " | Anneal:",
            start$par, " | ", e$message)
        if (!is.null(extra_info)) {
            cat(" | ", extra_info)
        }
        cat("\n")
    })
    theta.hat
}



exp_experiment_4_gen <- function(
    csv_filename,
    R = 1000,
    bernoulli_probs = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
    q = .25,
    sample_size = 100,
    append = TRUE,
    use_aneal_start = TRUE) {

    print(getwd())

    set.seed(32861) # set seed for reproducibility
    theta <- c(1,   # component 1 failure rate
            1.1,    # 2
            0.98,   # 3
            1.12,   # 4
            1.05)   # 5

    m <- length(theta)
    tau <- -(1/sum(theta))*log(q)

    if (!append) {
        cnames <- c("R", "p", "tau", "N", paste0("bias",1:m), "mse",
            paste0("se",1:m), paste0("se_asym",1:m), "mse_asym", "mse_asym_hat",
            paste0("coverage",1:m))

        # Write column names first
        write.table(t(cnames), file = csv_filename,
            sep = ",", col.names = FALSE,
            row.names = FALSE, append = FALSE, quote = FALSE)
    }

    for (i in 1:length(bernoulli_probs)) {
        p <- bernoulli_probs[i]
        cat("Starting simulations for Bernoulli probability", p, "\n")
        mles <- matrix(nrow = R, ncol = m)

        # For storing the lower and upper bounds of the confidence intervals
        CI_lwr <- matrix(nrow = R, ncol = m)
        CI_upr <- matrix(nrow = R, ncol = m)

        j <- 1L
        repeat {
            data <- generate_data(sample_size, theta, tau, p)
            theta.hat <- custom_solver(
                data = data,
                theta = theta,
                extra_info = paste0("Replicate(", j, ")"),
                annealing = use_aneal_start)
            if (is.null(theta.hat)) {
                next
            }
            if (j %% 10 == 0) {
                cat("Sample size", N, " | Replicate ", j,
                    " | MLE ", point(theta.hat), "\n")
            }
            mles[j, ] <- point(theta.hat)

            CI <- confint(theta.hat)
            if (any(is.nan(CI)))
            {
                print("NaN in CI")
                print(summary(theta.hat))
            }
            CI_lwr[j, ] <- CI[,1]
            CI_upr[j, ] <- CI[,2]

            j <- j + 1L
            if (j > R) {
                break
            }
        }

        # compute asymptotics
        theta.mle <- NULL
        while (is.null(theta.mle)) {
            data <- generate_data(sample_size, theta, tau, p)
            theta.mle <- custom_solver(
                data = data,
                theta = theta,
                extra_info = "asymptotics",
                annealing = use_aneal_start)
        }
        SE.asym <- se(theta.mle)
        MSE.asym <- mse(theta.hat, theta)
        MSE.asym.hat <- mse(theta.mle)

        # Calculate bias, MSE, and SE for each parameter
        bias <- colMeans(mles) - theta
        MSE <- sum(colMeans((mles - theta)^2))
        sigma <- cov(mles) * ((sample_size - 1) / sample_size)
        SE <- sqrt(diag(sigma))

        # Compute coverage probabilities
        coverage <- colMeans((CI_lwr <= theta) & (theta <= CI_upr))

        datum <- c(R, p, tau, N, bias, MSE, SE, SE.asym,
            MSE.asym, MSE.asym.hat, coverage)

        write.table(t(datum), file = csv_filename, sep = ",", col.names = FALSE,
            row.names = FALSE, append = TRUE)
    }
}


exp_experiment_4_gen(
    append = TRUE,
    csv_filename = "exp_experiment_4.csv",
    use_aneal_start = TRUE)

