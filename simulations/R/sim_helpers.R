# =============================================================================
# Simulation Helper Functions
# =============================================================================
# Shared utilities for the exponential series system simulation framework.
# Provides data generation, model fitting, checkpointing, progress display,
# and summary statistics.
#
# Dependencies: maskedcauses, md.tools, dplyr
#
# Usage: source("R/sim_helpers.R")
# =============================================================================

suppressPackageStartupMessages({
    library(maskedcauses)
    library(md.tools)
    library(dplyr)
})

# --- Model setup (created once, reused across all experiments) ---
MODEL  <- exp_series_md_c1_c2_c3()
SOLVER <- fit(MODEL)


# =============================================================================
# Data generation
# =============================================================================

#' Generate one masked series system dataset
#'
#' @param theta Numeric vector of true component failure rates
#' @param n     Integer sample size (number of systems)
#' @param q     Survival probability P(T > tau), determines censoring time
#' @param p     Masking probability for Bernoulli candidate sets
#' @return A data frame with columns t, delta, x1..xm (plus latent columns)
generate_data <- function(theta, n, q, p) {
    m <- length(theta)
    tau <- rep(-(1 / sum(theta)) * log(q), n)

    comp_times <- matrix(nrow = n, ncol = m)
    for (j in 1:m) comp_times[, j] <- rexp(n, theta[j])
    comp_times <- md_encode_matrix(comp_times, "t")

    comp_times %>%
        md_series_lifetime_right_censoring(tau) %>%
        md_bernoulli_cand_c1_c2_c3(p) %>%
        md_cand_sampler()
}


# =============================================================================
# Model fitting
# =============================================================================

#' Fit model to one dataset, returning structured results
#'
#' Returns a list with par, se, ci_lower, ci_upper, converged, loglik.
#' On failure, all numeric fields are NA and converged is FALSE.
#'
#' @param data   Data frame from generate_data()
#' @param theta0 Initial parameter guess (numeric vector of length m)
#' @param method Optimization method (default: "Nelder-Mead")
#' @param alpha  Significance level for Wald CIs (default: 0.05)
fit_one <- function(data, theta0, method = "Nelder-Mead", alpha = 0.05) {
    m <- length(theta0)
    na_result <- list(
        par       = rep(NA_real_, m),
        se        = rep(NA_real_, m),
        ci_lower  = rep(NA_real_, m),
        ci_upper  = rep(NA_real_, m),
        converged = FALSE,
        loglik    = NA_real_
    )

    tryCatch({
        result <- SOLVER(data, par = theta0, method = method)
        if (!result$converged) return(na_result)

        se <- sqrt(diag(result$vcov))
        z  <- qnorm(1 - alpha / 2)

        list(
            par       = result$par,
            se        = se,
            ci_lower  = result$par - z * se,
            ci_upper  = result$par + z * se,
            converged = TRUE,
            loglik    = result$loglik
        )
    }, error = function(e) na_result)
}


# =============================================================================
# Checkpointing (resumable simulations)
# =============================================================================

#' Load checkpoint from RDS file, or return NULL if none exists
load_checkpoint <- function(path) {
    if (file.exists(path)) readRDS(path) else NULL
}

#' Save checkpoint to RDS file (atomic write via tempfile + rename)
#'
#' Uses write-to-temp-then-rename to prevent corruption if the process
#' is interrupted during the write.
save_checkpoint <- function(data, path) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
    saveRDS(data, tmp)
    file.rename(tmp, path)
    invisible(path)
}


# =============================================================================
# Summary statistics
# =============================================================================

#' Compute MC summary statistics from an estimates matrix
#'
#' @param estimates Matrix (B x m) of parameter estimates, may contain NAs
#' @param theta     True parameter vector (length m)
#' @return Data frame with per-component bias, variance, MSE, RMSE, etc.
mc_summary <- function(estimates, theta) {
    m <- length(theta)
    valid <- complete.cases(estimates)
    n_valid <- sum(valid)

    if (n_valid < 2) {
        return(data.frame(
            component = 1:m, true = theta,
            mean_est = NA_real_, bias = NA_real_, variance = NA_real_,
            mse = NA_real_, rmse = NA_real_, rel_bias_pct = NA_real_,
            n_valid = n_valid
        ))
    }

    est <- estimates[valid, , drop = FALSE]
    bias     <- colMeans(est) - theta
    variance <- apply(est, 2, var)
    mse      <- bias^2 + variance

    data.frame(
        component    = 1:m,
        true         = theta,
        mean_est     = colMeans(est),
        bias         = bias,
        variance     = variance,
        mse          = mse,
        rmse         = sqrt(mse),
        rel_bias_pct = 100 * bias / theta,
        n_valid      = n_valid
    )
}

#' Compute CI coverage statistics
#'
#' @param ci_lower Matrix (B x m) of CI lower bounds
#' @param ci_upper Matrix (B x m) of CI upper bounds
#' @param theta    True parameter vector
#' @param alpha    Significance level
#' @return Data frame with per-component coverage and mean CI width
coverage_summary <- function(ci_lower, ci_upper, theta, alpha = 0.05) {
    m <- length(theta)
    valid <- complete.cases(ci_lower) & complete.cases(ci_upper)
    n_valid <- sum(valid)

    if (n_valid < 1) {
        return(data.frame(
            component = 1:m, true = theta,
            coverage = NA_real_, nominal = 1 - alpha,
            mean_width = NA_real_, n_valid = 0L
        ))
    }

    coverage   <- numeric(m)
    mean_width <- numeric(m)
    for (j in 1:m) {
        covered        <- ci_lower[valid, j] <= theta[j] &
                          theta[j] <= ci_upper[valid, j]
        coverage[j]    <- mean(covered)
        mean_width[j]  <- mean(ci_upper[valid, j] - ci_lower[valid, j])
    }

    data.frame(
        component  = 1:m,
        true       = theta,
        coverage   = coverage,
        nominal    = 1 - alpha,
        mean_width = mean_width,
        n_valid    = n_valid
    )
}

#' Compute sensitivity summary (per-component + means across components)
#'
#' @param estimates Matrix (B x m) of estimates, may contain NAs
#' @param theta     True parameter vector
#' @return List with per-component and mean summary statistics
sensitivity_summary <- function(estimates, theta) {
    valid <- complete.cases(estimates)
    n_valid <- sum(valid)

    if (n_valid < 2) {
        return(list(
            bias = rep(NA_real_, length(theta)),
            variance = rep(NA_real_, length(theta)),
            mse = rep(NA_real_, length(theta)),
            mean_abs_bias = NA_real_,
            mean_mse = NA_real_,
            mean_rmse = NA_real_,
            n_valid = n_valid
        ))
    }

    est      <- estimates[valid, , drop = FALSE]
    bias     <- colMeans(est) - theta
    variance <- apply(est, 2, var)
    mse      <- bias^2 + variance

    list(
        bias          = bias,
        variance      = variance,
        mse           = mse,
        mean_abs_bias = mean(abs(bias)),
        mean_mse      = mean(mse),
        mean_rmse     = sqrt(mean(mse)),
        n_valid       = n_valid
    )
}


# =============================================================================
# Progress display
# =============================================================================

#' Print a timestamped log message
log_msg <- function(...) {
    cat(format(Sys.time(), "[%H:%M:%S] "), paste0(...), "\n", sep = "")
    flush.console()
}

#' Display progress for a replication loop
#'
#' Call this at checkpoint intervals (every SAVE_EVERY reps).
#' Shows: [current/total] percent | rate | ETA
#'
#' @param b     Current replication number
#' @param total Total replications
#' @param t0    Start time from proc.time()["elapsed"]
#' @param start First replication in this run (for rate calculation)
#' @param label Optional prefix label
log_progress <- function(b, total, t0, start = 1L, label = "") {
    elapsed <- proc.time()["elapsed"] - t0
    done    <- b - start + 1L
    rate    <- if (elapsed > 0) done / elapsed else 0
    remaining <- if (rate > 0) (total - b) / rate else Inf
    pct     <- 100 * b / total

    prefix <- if (nzchar(label)) paste0("  ", label, " ") else "  "
    cat(sprintf("%s[%*d/%d] %5.1f%% | %.1f reps/s | ETA: %s\n",
                prefix, nchar(as.character(total)), b, total, pct,
                rate, format_elapsed(remaining)))
    flush.console()
}

#' Format elapsed seconds as human-readable string
format_elapsed <- function(secs) {
    if (is.na(secs) || !is.finite(secs)) return("?")
    if (secs < 60) return(sprintf("%.0fs", secs))
    if (secs < 3600) return(sprintf("%.1fmin", secs / 60))
    sprintf("%.1fh", secs / 3600)
}

#' Check for --reset flag in command line arguments
should_reset <- function() {
    "--reset" %in% commandArgs(trailingOnly = TRUE)
}
