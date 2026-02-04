# =============================================================================
# Experiment 4: Sample Size Study (with Coverage Analysis)
# =============================================================================
# Studies how sample size n affects MLE properties: bias, variance, MSE,
# CI coverage, and CI width. Verifies the 1/sqrt(n) scaling of CI width
# and shows at what n the asymptotic Wald approximation becomes reliable.
#
# Parameters: n in {100, 250, 500, 1000, 2500, 5000}, B=100, p=0.3, ~25% cens
# Resumable:  saves checkpoint every SAVE_EVERY replications
#
# Usage:
#   cd simulations
#   Rscript experiments/04_sample_size.R
# =============================================================================

source("R/sim_helpers.R")
source("experiments/config.R")

RESULT_FILE <- "results/04_sample_size.rds"
n_sizes <- length(SIZE_N_VALUES)

if (should_reset() && file.exists(RESULT_FILE)) {
    file.remove(RESULT_FILE)
    log_msg("Reset: removed ", RESULT_FILE)
}

# --- Load or initialize checkpoint ---
cp <- load_checkpoint(RESULT_FILE)
if (is.null(cp)) {
    cp <- list(
        config       = list(n_values = SIZE_N_VALUES, B = SIZE_B,
                            p = SIZE_P, q = SIZE_Q,
                            theta = THETA, alpha = ALPHA, seed = SEED),
        estimates    = lapply(1:n_sizes, function(i) matrix(NA_real_, SIZE_B, M)),
        se_estimates = lapply(1:n_sizes, function(i) matrix(NA_real_, SIZE_B, M)),
        ci_lower     = lapply(1:n_sizes, function(i) matrix(NA_real_, SIZE_B, M)),
        ci_upper     = lapply(1:n_sizes, function(i) matrix(NA_real_, SIZE_B, M)),
        converged    = lapply(1:n_sizes, function(i) rep(NA, SIZE_B)),
        completed    = integer(n_sizes)
    )
}

# --- Run ---
total_fits <- sum(pmax(0, SIZE_B - cp$completed))
log_msg("Experiment 4: Sample size study")
log_msg("  n in {", paste(SIZE_N_VALUES, collapse = ", "), "}")
log_msg("  B=", SIZE_B, ", p=", SIZE_P, ", q=", SIZE_Q)
log_msg("  ", total_fits, " fits remaining")

t0_total <- proc.time()["elapsed"]

for (n_idx in seq_along(SIZE_N_VALUES)) {
    n_curr  <- SIZE_N_VALUES[n_idx]
    start_b <- cp$completed[n_idx] + 1L

    if (start_b > SIZE_B) next

    log_msg("  n=", n_curr,
            if (start_b > 1L) paste0(" (resuming from rep ", start_b, ")")
            else "")
    t0 <- proc.time()["elapsed"]

    for (b in start_b:SIZE_B) {
        set.seed(SEED + 300000L + n_idx * 1000L + b)

        data_b <- generate_data(THETA, n_curr, SIZE_Q, SIZE_P)
        result <- fit_one(data_b, THETA0, METHOD, ALPHA)

        cp$estimates[[n_idx]][b, ]    <- result$par
        cp$se_estimates[[n_idx]][b, ] <- result$se
        cp$ci_lower[[n_idx]][b, ]     <- result$ci_lower
        cp$ci_upper[[n_idx]][b, ]     <- result$ci_upper
        cp$converged[[n_idx]][b]      <- result$converged
        cp$completed[n_idx]           <- b

        if (b %% SAVE_EVERY == 0 || b == SIZE_B) {
            save_checkpoint(cp, RESULT_FILE)
            log_progress(b, SIZE_B, t0, start_b,
                         label = sprintf("n=%d", n_curr))
        }
    }
}

total_elapsed <- proc.time()["elapsed"] - t0_total
log_msg("Experiment 4 complete in ", format_elapsed(total_elapsed))

# --- Print summary ---
cat("\n")
log_msg("=== Sample Size Study Summary ===")

size_df <- data.frame(
    n         = SIZE_N_VALUES,
    conv_rate = sapply(seq_along(SIZE_N_VALUES), function(i)
        mean(cp$converged[[i]], na.rm = TRUE)),
    mean_rmse = sapply(seq_along(SIZE_N_VALUES), function(i)
        mc_summary(cp$estimates[[i]], THETA)$rmse |> mean()),
    mean_coverage = sapply(seq_along(SIZE_N_VALUES), function(i)
        coverage_summary(cp$ci_lower[[i]], cp$ci_upper[[i]],
                         THETA, ALPHA)$coverage |> mean()),
    mean_ci_width = sapply(seq_along(SIZE_N_VALUES), function(i)
        coverage_summary(cp$ci_lower[[i]], cp$ci_upper[[i]],
                         THETA, ALPHA)$mean_width |> mean())
)
print(size_df, digits = 4, row.names = FALSE)

# Verify 1/sqrt(n) scaling
cat("\nCI width ratio vs 1/sqrt(n) ratio (should be ~1.0 if scaling holds):\n")
ref_idx <- which(SIZE_N_VALUES == 500)
if (length(ref_idx) == 1) {
    for (i in seq_along(SIZE_N_VALUES)) {
        width_ratio <- size_df$mean_ci_width[i] / size_df$mean_ci_width[ref_idx]
        theory_ratio <- sqrt(SIZE_N_VALUES[ref_idx] / SIZE_N_VALUES[i])
        cat(sprintf("  n=%5d: width_ratio=%.3f, 1/sqrt(n) ratio=%.3f, ratio=%.3f\n",
                    SIZE_N_VALUES[i], width_ratio, theory_ratio,
                    width_ratio / theory_ratio))
    }
}
