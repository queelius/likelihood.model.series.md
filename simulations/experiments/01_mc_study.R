# =============================================================================
# Experiment 1: Main Monte Carlo Simulation Study
# =============================================================================
# Estimates bias, variance, MSE, RMSE, CI coverage, and CI width of the MLE
# for the exponential series system with masked data.
#
# Parameters: n=7500, B=200, p=0.3, ~25% censoring
# Resumable:  saves checkpoint every SAVE_EVERY replications
# Reset:      Rscript experiments/01_mc_study.R --reset
#
# Usage:
#   cd simulations
#   Rscript experiments/01_mc_study.R
# =============================================================================

source("R/sim_helpers.R")
source("experiments/config.R")

RESULT_FILE <- "results/01_mc_study.rds"

# --- Handle --reset flag ---
if (should_reset() && file.exists(RESULT_FILE)) {
    file.remove(RESULT_FILE)
    log_msg("Reset: removed ", RESULT_FILE)
}

# --- Load or initialize checkpoint ---
cp <- load_checkpoint(RESULT_FILE)
if (is.null(cp)) {
    cp <- list(
        config       = list(n = MC_N, B = MC_B, p = MC_P, q = MC_Q,
                            theta = THETA, alpha = ALPHA, seed = SEED),
        estimates    = matrix(NA_real_, MC_B, M),
        se_estimates = matrix(NA_real_, MC_B, M),
        ci_lower     = matrix(NA_real_, MC_B, M),
        ci_upper     = matrix(NA_real_, MC_B, M),
        converged    = rep(NA, MC_B),
        completed    = 0L
    )
}

# --- Run remaining replications ---
start_b <- cp$completed + 1L

if (start_b > MC_B) {
    log_msg("Experiment 1 already complete (", MC_B, "/", MC_B, " reps)")
} else {
    log_msg("Experiment 1: MC study")
    log_msg("  n=", MC_N, ", B=", MC_B, ", p=", MC_P,
            ", q=", MC_Q, " (~", round(100 * MC_Q), "% censoring)")

    if (start_b > 1L) {
        log_msg("  Resuming from rep ", start_b, "/", MC_B,
                " (", start_b - 1L, " already done)")
    }

    t0 <- proc.time()["elapsed"]

    for (b in start_b:MC_B) {
        set.seed(SEED + b)

        data_b <- generate_data(THETA, MC_N, MC_Q, MC_P)
        result <- fit_one(data_b, THETA0, METHOD, ALPHA)

        cp$estimates[b, ]    <- result$par
        cp$se_estimates[b, ] <- result$se
        cp$ci_lower[b, ]     <- result$ci_lower
        cp$ci_upper[b, ]     <- result$ci_upper
        cp$converged[b]      <- result$converged
        cp$completed         <- b

        if (b %% SAVE_EVERY == 0 || b == MC_B) {
            save_checkpoint(cp, RESULT_FILE)
            log_progress(b, MC_B, t0, start_b)
        }
    }

    elapsed <- proc.time()["elapsed"] - t0
    log_msg("Experiment 1 complete in ", format_elapsed(elapsed))
}

# --- Print summary ---
cat("\n")
log_msg("=== MC Study Summary ===")
cat("\nConvergence rate:", mean(cp$converged, na.rm = TRUE), "\n")

cat("\nBias, Variance, MSE:\n")
print(mc_summary(cp$estimates, THETA), digits = 4, row.names = FALSE)

cat("\nCI Coverage:\n")
print(coverage_summary(cp$ci_lower, cp$ci_upper, THETA, ALPHA),
      digits = 4, row.names = FALSE)
