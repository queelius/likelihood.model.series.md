# =============================================================================
# Experiment 5: Joint Masking x Censoring Interaction
# =============================================================================
# Studies whether masking and censoring effects on MLE accuracy are additive
# or interact nonlinearly. Runs a grid of (p, q) combinations.
#
# Parameters: n=500, B=100, p in {0, 0.2, 0.4}, q in {0.1, 0.3, 0.5, 0.7, 0.9}
# Grid:       3 x 5 = 15 cells
# Resumable:  saves checkpoint every SAVE_EVERY replications
#
# Usage:
#   cd simulations
#   Rscript experiments/05_joint_masking_censoring.R
# =============================================================================

source("R/sim_helpers.R")
source("experiments/config.R")

RESULT_FILE <- "results/05_joint_masking_censoring.rds"
n_p <- length(JOINT_P_VALUES)
n_q <- length(JOINT_Q_VALUES)

if (should_reset() && file.exists(RESULT_FILE)) {
    file.remove(RESULT_FILE)
    log_msg("Reset: removed ", RESULT_FILE)
}

# --- Load or initialize checkpoint ---
cp <- load_checkpoint(RESULT_FILE)
if (is.null(cp)) {
    cp <- list(
        config     = list(n = JOINT_N, B = JOINT_B,
                          p_values = JOINT_P_VALUES,
                          q_values = JOINT_Q_VALUES,
                          theta = THETA, seed = SEED),
        # Nested list: estimates[[p_idx]][[q_idx]] is B x M matrix
        estimates  = lapply(1:n_p, function(i)
            lapply(1:n_q, function(j) matrix(NA_real_, JOINT_B, M))),
        cens_rates = lapply(1:n_p, function(i)
            lapply(1:n_q, function(j) rep(NA_real_, JOINT_B))),
        completed  = matrix(0L, n_p, n_q)
    )
}

# --- Run ---
total_fits <- sum(pmax(0L, JOINT_B - cp$completed))
log_msg("Experiment 5: Joint masking x censoring")
log_msg("  n=", JOINT_N, ", B=", JOINT_B)
log_msg("  p in {", paste(JOINT_P_VALUES, collapse = ", "), "}")
log_msg("  q in {", paste(JOINT_Q_VALUES, collapse = ", "), "}")
log_msg("  Grid: ", n_p, " x ", n_q, " = ", n_p * n_q, " cells")
log_msg("  ", total_fits, " fits remaining")

t0_total <- proc.time()["elapsed"]
cell_num <- 0L
total_cells <- n_p * n_q

for (p_idx in seq_len(n_p)) {
    p_curr <- JOINT_P_VALUES[p_idx]

    for (q_idx in seq_len(n_q)) {
        q_curr  <- JOINT_Q_VALUES[q_idx]
        start_b <- cp$completed[p_idx, q_idx] + 1L
        cell_num <- cell_num + 1L

        if (start_b > JOINT_B) next

        log_msg("  Cell ", cell_num, "/", total_cells,
                ": p=", sprintf("%.1f", p_curr),
                ", q=", sprintf("%.1f", q_curr),
                if (start_b > 1L) paste0(" (resuming from rep ", start_b, ")")
                else "")
        t0 <- proc.time()["elapsed"]

        for (b in start_b:JOINT_B) {
            set.seed(SEED + 400000L + p_idx * 100000L +
                     q_idx * 1000L + b)

            data_b <- generate_data(THETA, JOINT_N, q_curr, p_curr)
            result <- fit_one(data_b, THETA0, METHOD, ALPHA)

            if (result$converged) {
                cp$estimates[[p_idx]][[q_idx]][b, ] <- result$par
            }
            cp$cens_rates[[p_idx]][[q_idx]][b] <- 1 - mean(data_b$delta)
            cp$completed[p_idx, q_idx] <- b

            if (b %% SAVE_EVERY == 0 || b == JOINT_B) {
                save_checkpoint(cp, RESULT_FILE)
                log_progress(b, JOINT_B, t0, start_b,
                             label = sprintf("p=%.1f,q=%.1f", p_curr, q_curr))
            }
        }
    }
}

total_elapsed <- proc.time()["elapsed"] - t0_total
log_msg("Experiment 5 complete in ", format_elapsed(total_elapsed))

# --- Print summary (MSE heat map as table) ---
cat("\n")
log_msg("=== Joint Masking x Censoring: Mean MSE ===")
cat("\n")

# Column headers
cat(sprintf("%8s", "p \\ q"))
for (q_idx in seq_len(n_q)) {
    cat(sprintf(" %7.1f", JOINT_Q_VALUES[q_idx]))
}
cat("\n")

for (p_idx in seq_len(n_p)) {
    cat(sprintf("%8.1f", JOINT_P_VALUES[p_idx]))
    for (q_idx in seq_len(n_q)) {
        ss <- sensitivity_summary(cp$estimates[[p_idx]][[q_idx]], THETA)
        cat(sprintf(" %7.4f", ss$mean_mse))
    }
    cat("\n")
}
