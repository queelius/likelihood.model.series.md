library(utils)
library(ggplot2)
library(tidyverse)
library(algebraic.mle)

df <- read.csv("exp_experiment_3.csv")

###################### coverage probability ##################
df_long <- df %>% pivot_longer(
    cols = starts_with("coverage"),
    names_to = "Coverage",
    values_to = "Value")

# Convert the "Coverage" column to a factor for better plotting
df_long$Coverage <- as.factor(df_long$Coverage)

# Plot
ggplot(df_long, aes(x = N, y = Value, color = Coverage)) +
  geom_point() +
  geom_line() +
  labs(x = "Sample Size (N)", y = "Coverage Probability",
    title = "Coverage Probability vs Sample Size", color = "Coverage")
  



############# BIAS #############
df_long <- df %>% pivot_longer(
    cols = starts_with("bias"),
    names_to = "Bias",
    values_to = "Value")

# Convert the "Coverage" column to a factor for better plotting
df_long$Bias <- as.factor(df_long$Bias)

# Plot
ggplot(df_long, aes(x = N, y = Value, color = Bias)) +
  geom_point() +
  geom_line() +
  labs(x = "Sample Size (N)", y = "Bias", title = "Bias vs Sample Size", color = "Bias")
  #theme_minimal()


ggplot(df) +
  geom_line(aes(x=N, y=bias1, color="Bias(1)")) +
  geom_line(aes(x=N, y=bias2, color="Bias(2)")) +
  geom_line(aes(x=N, y=bias3, color="Bias(3)")) +
  geom_line(aes(x=N, y=bias4, color="Bias(4)")) +
  geom_line(aes(x=N, y=bias5, color="Bias(5)")) +
  geom_line(aes(x=N, y=bias6, color="Bias(6)")) +
  geom_line(aes(x=N, y=bias7, color="Bias(7)")) +
  labs(x = "Sample Size (N)", y = "Bias", title = "Bias vs Sample Size")

############### MSE ###############

ggplot(df) +
  geom_point(aes(x=N, y=mse)) +
  geom_line(aes(x=N, y=mse, color = "Simulation")) +
  geom_line(aes(x=N,y=mse_asym,color="Asymptotic"),linetype="dashed") +
  geom_line(aes(x=N,y=mse_asym_hat,color="Estimated Asymptotic"),linetype="dashed") +
  labs(x = "Sample Size (N)", y = "Mean Squared Error", title = "Mean Squared Error vs Sample Size")


############## SE ###############
df_long <- df %>% pivot_longer(
    # regex to match column names "se<digit>"
    cols = matches("se[0-9]"),
    names_to = "SE",
    values_to = "Value")

# Convert the "Coverage" column to a factor for better plotting
df_long$SE <- as.factor(df_long$SE)

# Plot
ggplot(df_long, aes(x = N, y = Value, color = SE)) +
  geom_point() +
  geom_line() +
  labs(x = "Sample Size (N)", y = "SE", title = "SE vs Sample Size", color = "SE")
  #theme_minimal()

############## SE asymptotics ###############
df_long <- df %>% pivot_longer(
    # regex to match column names "se<digit>"
    cols = starts_with("se_asym"),
    names_to = "SE",
    values_to = "Value")

# Convert the "Coverage" column to a factor for better plotting
df_long$SE <- as.factor(df_long$SE)

# Plot
ggplot(df_long, aes(x = N, y = Value, color = SE)) +
  geom_point() +
  geom_line() +
  labs(x = "Sample Size (N)", y = "SE", title = "SE vs Sample Size", color = "SE")
  #theme_minimal()


#############################################################################

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

exp_experiment_3_gen <- function(
    csv_filename,
    R = 1000,
    p = .333,
    q = .25,
    sample_sizes = c(
        5, 10, 15, 20, 25, 30, 35, 40,
        45, 50, 60, 70, 80, 90, 100, 250, 500),
    append = TRUE,
    use_aneal_start = TRUE) {
    # 7-component series system
    # let's think about failure rates. with the weibull, we can make it so
    # that some components fail early or last a long time,
    # others the reverse, and the remaining somewhere inbetween, so we get
    # a nice distribution of component failures.


    # TODO: get one or more of the samples that fail to converge (gradient
    # fn doens't work), and investigate it in more detail. you my find
    # out something interesting about the data that way.

    set.seed(7231) # set seed for reproducibility
    theta <- c(1,  # component 1 failure rate
            1.1,   # 2
            0.975, # 3
            1.125, # 4
            1.1,   # 5
            1.0,   # 6
            1.05)  # 7

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

    for (i in 1:length(sample_sizes)) {
        N <- sample_sizes[i]
        cat("Starting simulations for sample size", N, "\n")
        mles <- matrix(nrow = R, ncol = m)

        # For storing the lower and upper bounds of the confidence intervals
        CI_lwr <- matrix(nrow = R, ncol = m)
        CI_upr <- matrix(nrow = R, ncol = m)

        j <- 1L
        repeat {
            data <- generate_data(N, theta, tau, p)
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
            data <- generate_data(N, theta, tau, p)
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
        sigma <- cov(mles) * ((N - 1) / N) # MLE of variance-covariance matrix
        SE <- sqrt(diag(sigma))

        # Compute coverage probabilities
        coverage <- colMeans((CI_lwr <= theta) & (theta <= CI_upr))

        datum <- c(R, p, tau, N, bias, MSE, SE, SE.asym,
            MSE.asym, MSE.asym.hat, coverage)

        write.table(t(datum), file = csv_filename, sep = ",", col.names = FALSE,
            row.names = FALSE, append = TRUE)
    }
}

exp_experiment_3_gen(
    csv_filename = "exp_experiment_3-1.csv",
    sample_sizes = c(1000, 2000),
    #R = 1000, p = .333, q = .25,
    append = FALSE,
    use_aneal_start = TRUE)

###########################
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

