# =============================================================================
# Export Simulation Results for Vignette
# =============================================================================
# Transforms the simulation framework's results into the format expected
# by vignettes/exponential_series.Rmd (which loads via list2env).
#
# The vignette expects these variables:
#   estimates, se_estimates, ci_lower, ci_upper, converged  (MC study)
#   B, alpha                                                (MC config)
#   mask_results, p_values                                  (masking sens.)
#   cens_results, q_values                                  (censoring sens.)
#   theta, m, n_sens, B_sens, p_cens                        (shared config)
#
# Additional results (sample size, joint) are included under separate
# keys for future vignette expansion.
#
# Usage:
#   cd simulations
#   Rscript export_for_vignette.R
# =============================================================================

source("R/sim_helpers.R")
source("experiments/config.R")

# --- Check that required results exist ---
required <- c("results/01_mc_study.rds",
              "results/02_masking_sensitivity.rds",
              "results/03_censoring_sensitivity.rds")
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
    stop("Missing required results:\n  ",
         paste(missing, collapse = "\n  "),
         "\nRun the experiments first: Rscript run_all.R")
}

# --- Load results ---
mc   <- readRDS("results/01_mc_study.rds")
mask <- readRDS("results/02_masking_sensitivity.rds")
cens <- readRDS("results/03_censoring_sensitivity.rds")

# --- Transform MC study (Experiment 1) ---
# These become top-level variables in the vignette environment
estimates    <- mc$estimates
se_estimates <- mc$se_estimates
ci_lower     <- mc$ci_lower
ci_upper     <- mc$ci_upper
converged    <- mc$converged
B            <- mc$config$B
alpha_val    <- mc$config$alpha

# --- Transform masking sensitivity (Experiment 2) ---
# Format: list of lists, each with p, bias, variance, mse vectors
mask_results <- lapply(seq_along(mask$config$p_values), function(i) {
    ss <- sensitivity_summary(mask$estimates[[i]], mask$config$theta)
    list(
        p        = mask$config$p_values[i],
        bias     = ss$bias,
        variance = ss$variance,
        mse      = ss$mse
    )
})
p_values <- mask$config$p_values

# --- Transform censoring sensitivity (Experiment 3) ---
# Format: list of lists, each with q, cens_rate, bias, variance, mse
cens_results <- lapply(seq_along(cens$config$q_values), function(i) {
    ss <- sensitivity_summary(cens$estimates[[i]], cens$config$theta)
    list(
        q         = cens$config$q_values[i],
        cens_rate = mean(cens$cens_rates[[i]], na.rm = TRUE),
        bias      = ss$bias,
        variance  = ss$variance,
        mse       = ss$mse
    )
})
q_values <- cens$config$q_values

# --- Assemble the export ---
export <- list(
    # MC study
    estimates    = estimates,
    se_estimates = se_estimates,
    ci_lower     = ci_lower,
    ci_upper     = ci_upper,
    converged    = converged,
    B            = B,
    alpha        = alpha_val,

    # Masking sensitivity
    mask_results = mask_results,
    p_values     = p_values,

    # Censoring sensitivity
    cens_results = cens_results,
    q_values     = q_values,

    # Shared config
    theta  = mc$config$theta,
    m      = length(mc$config$theta),
    n_sens = mask$config$n,
    B_sens = mask$config$B,
    p_cens = cens$config$p
)

# --- Add optional experiment results ---
if (file.exists("results/04_sample_size.rds")) {
    export$sample_size_results <- readRDS("results/04_sample_size.rds")
    log_msg("  Included: sample size study")
}

if (file.exists("results/05_joint_masking_censoring.rds")) {
    export$joint_results <- readRDS("results/05_joint_masking_censoring.rds")
    log_msg("  Included: joint masking x censoring study")
}

# --- Save ---
out_path <- file.path("..", "vignettes", "precomputed_results.rds")
saveRDS(export, out_path)
log_msg("Exported to: ", normalizePath(out_path))

# --- Verification ---
cat("\nExport verification:\n")
cat("  MC study:   ", nrow(estimates), " reps,",
    sum(converged, na.rm = TRUE), "converged\n")
cat("  Masking:    ", length(p_values), " p values\n")
cat("  Censoring:  ", length(q_values), " q values\n")
cat("  theta:      ", paste(export$theta, collapse = ", "), "\n")
cat("  m:          ", export$m, "\n")
