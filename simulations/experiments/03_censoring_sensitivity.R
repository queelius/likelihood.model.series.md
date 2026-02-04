# =============================================================================
# Experiment 3: Right-Censoring Sensitivity
# =============================================================================
# Studies how the right-censoring rate affects MLE accuracy.
# Varies the survival probability q from 0.1 (light censoring) to 0.9
# (heavy censoring). Fixed masking at p=0.2.
#
# Improvement over original: 9 q values (vs 5) for smoother curves.
# Also tracks actual empirical censoring rate per replication.
#
# Parameters: n=500, B=100 per q value, q in {0.1, 0.2, ..., 0.9}
# Resumable:  saves checkpoint every SAVE_EVERY replications
#
# Usage:
#   cd simulations
#   Rscript experiments/03_censoring_sensitivity.R
# =============================================================================

source("R/sim_helpers.R")
source("experiments/config.R")

RESULT_FILE <- "results/03_censoring_sensitivity.rds"
n_q <- length(CENS_Q_VALUES)

if (should_reset() && file.exists(RESULT_FILE)) {
    file.remove(RESULT_FILE)
    log_msg("Reset: removed ", RESULT_FILE)
}

# --- Load or initialize checkpoint ---
cp <- load_checkpoint(RESULT_FILE)
if (is.null(cp)) {
    cp <- list(
        config     = list(n = SENS_N, B = SENS_B, q_values = CENS_Q_VALUES,
                          p = CENS_P, theta = THETA, seed = SEED),
        estimates  = lapply(1:n_q, function(i) matrix(NA_real_, SENS_B, M)),
        cens_rates = lapply(1:n_q, function(i) rep(NA_real_, SENS_B)),
        completed  = integer(n_q)
    )
}

# --- Run ---
total_fits <- sum(pmax(0, SENS_B - cp$completed))
log_msg("Experiment 3: Censoring sensitivity")
log_msg("  n=", SENS_N, ", B=", SENS_B, ", p=", CENS_P,
        ", q in {", paste(CENS_Q_VALUES, collapse = ", "), "}")
log_msg("  ", total_fits, " fits remaining")

t0_total <- proc.time()["elapsed"]

for (q_idx in seq_along(CENS_Q_VALUES)) {
    q_curr  <- CENS_Q_VALUES[q_idx]
    start_b <- cp$completed[q_idx] + 1L

    if (start_b > SENS_B) next

    log_msg("  q=", sprintf("%.1f", q_curr),
            " (theoretical ~", round(100 * q_curr), "% censoring)",
            if (start_b > 1L) paste0(" resuming from rep ", start_b)
            else "")
    t0 <- proc.time()["elapsed"]

    for (b in start_b:SENS_B) {
        set.seed(SEED + 200000L + q_idx * 1000L + b)

        data_b <- generate_data(THETA, SENS_N, q_curr, CENS_P)
        result <- fit_one(data_b, THETA0, METHOD, ALPHA)

        if (result$converged) {
            cp$estimates[[q_idx]][b, ] <- result$par
        }
        cp$cens_rates[[q_idx]][b] <- 1 - mean(data_b$delta)
        cp$completed[q_idx] <- b

        if (b %% SAVE_EVERY == 0 || b == SENS_B) {
            save_checkpoint(cp, RESULT_FILE)
            log_progress(b, SENS_B, t0, start_b,
                         label = sprintf("q=%.1f", q_curr))
        }
    }
}

total_elapsed <- proc.time()["elapsed"] - t0_total
log_msg("Experiment 3 complete in ", format_elapsed(total_elapsed))

# --- Print summary ---
cat("\n")
log_msg("=== Censoring Sensitivity Summary ===")
cens_df <- data.frame(
    q         = CENS_Q_VALUES,
    cens_pct  = sapply(seq_along(CENS_Q_VALUES), function(i)
        round(100 * mean(cp$cens_rates[[i]], na.rm = TRUE), 1)),
    mean_abs_bias = sapply(seq_along(CENS_Q_VALUES), function(i)
        sensitivity_summary(cp$estimates[[i]], THETA)$mean_abs_bias),
    mean_mse      = sapply(seq_along(CENS_Q_VALUES), function(i)
        sensitivity_summary(cp$estimates[[i]], THETA)$mean_mse),
    mean_rmse     = sapply(seq_along(CENS_Q_VALUES), function(i)
        sensitivity_summary(cp$estimates[[i]], THETA)$mean_rmse)
)
print(cens_df, digits = 4, row.names = FALSE)
