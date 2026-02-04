# =============================================================================
# Experiment 2: Masking Probability Sensitivity
# =============================================================================
# Studies how the masking probability p affects MLE accuracy.
# Varies p from 0 (perfect component ID) to 0.5 (maximum Bernoulli masking).
# Fixed censoring at ~25%.
#
# Parameters: n=500, B=100 per p value, p in {0, 0.1, ..., 0.5}
# Resumable:  saves checkpoint every SAVE_EVERY replications
#
# Usage:
#   cd simulations
#   Rscript experiments/02_masking_sensitivity.R
# =============================================================================

source("R/sim_helpers.R")
source("experiments/config.R")

RESULT_FILE <- "results/02_masking_sensitivity.rds"
n_p <- length(MASK_P_VALUES)

if (should_reset() && file.exists(RESULT_FILE)) {
    file.remove(RESULT_FILE)
    log_msg("Reset: removed ", RESULT_FILE)
}

# --- Load or initialize checkpoint ---
cp <- load_checkpoint(RESULT_FILE)
if (is.null(cp)) {
    cp <- list(
        config    = list(n = SENS_N, B = SENS_B, p_values = MASK_P_VALUES,
                         q = MASK_Q, theta = THETA, seed = SEED),
        estimates = lapply(1:n_p, function(i) matrix(NA_real_, SENS_B, M)),
        completed = integer(n_p)
    )
}

# --- Run ---
total_fits <- sum(pmax(0, SENS_B - cp$completed))
log_msg("Experiment 2: Masking sensitivity")
log_msg("  n=", SENS_N, ", B=", SENS_B, ", p in {",
        paste(MASK_P_VALUES, collapse = ", "), "}")
log_msg("  ", total_fits, " fits remaining")

t0_total <- proc.time()["elapsed"]

for (p_idx in seq_along(MASK_P_VALUES)) {
    p_curr  <- MASK_P_VALUES[p_idx]
    start_b <- cp$completed[p_idx] + 1L

    if (start_b > SENS_B) next

    log_msg("  p=", sprintf("%.1f", p_curr),
            if (start_b > 1L) paste0(" (resuming from rep ", start_b, ")")
            else "")
    t0 <- proc.time()["elapsed"]

    for (b in start_b:SENS_B) {
        set.seed(SEED + 100000L + p_idx * 1000L + b)

        data_b <- generate_data(THETA, SENS_N, MASK_Q, p_curr)
        result <- fit_one(data_b, THETA0, METHOD, ALPHA)

        if (result$converged) {
            cp$estimates[[p_idx]][b, ] <- result$par
        }
        cp$completed[p_idx] <- b

        if (b %% SAVE_EVERY == 0 || b == SENS_B) {
            save_checkpoint(cp, RESULT_FILE)
            log_progress(b, SENS_B, t0, start_b,
                         label = sprintf("p=%.1f", p_curr))
        }
    }
}

total_elapsed <- proc.time()["elapsed"] - t0_total
log_msg("Experiment 2 complete in ", format_elapsed(total_elapsed))

# --- Print summary ---
cat("\n")
log_msg("=== Masking Sensitivity Summary ===")
mask_df <- data.frame(
    p             = MASK_P_VALUES,
    mean_abs_bias = sapply(seq_along(MASK_P_VALUES), function(i)
        sensitivity_summary(cp$estimates[[i]], THETA)$mean_abs_bias),
    mean_mse      = sapply(seq_along(MASK_P_VALUES), function(i)
        sensitivity_summary(cp$estimates[[i]], THETA)$mean_mse),
    mean_rmse     = sapply(seq_along(MASK_P_VALUES), function(i)
        sensitivity_summary(cp$estimates[[i]], THETA)$mean_rmse)
)
print(mask_df, digits = 4, row.names = FALSE)
