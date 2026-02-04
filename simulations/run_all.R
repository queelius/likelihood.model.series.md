# =============================================================================
# Run All Experiments
# =============================================================================
# Orchestrator for the exponential series system simulation suite.
# Runs all experiments sequentially, then exports results for the vignette.
#
# Usage:
#   cd simulations
#   Rscript run_all.R              # Run all experiments
#   Rscript run_all.R 1 3 5        # Run only experiments 1, 3, and 5
#   Rscript run_all.R --reset      # Reset all experiments and re-run
#   Rscript run_all.R 1 --reset    # Reset and re-run experiment 1 only
#
# Each experiment is resumable: if interrupted, re-running picks up where
# it left off (unless --reset is passed).
# =============================================================================

# Verify working directory
if (!file.exists("R/sim_helpers.R")) {
    stop("Run this script from the simulations/ directory:\n",
         "  cd simulations && Rscript run_all.R\n",
         "Or from project root:\n",
         "  Rscript simulations/run_all.R")
}

experiments <- c(
    "experiments/01_mc_study.R",
    "experiments/02_masking_sensitivity.R",
    "experiments/03_censoring_sensitivity.R",
    "experiments/04_sample_size.R",
    "experiments/05_joint_masking_censoring.R"
)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
numeric_args <- suppressWarnings(as.integer(args[!args %in% "--reset"]))
numeric_args <- numeric_args[!is.na(numeric_args)]

if (length(numeric_args) > 0) {
    experiments <- experiments[numeric_args]
    cat("Running selected experiments:", paste(numeric_args, collapse = ", "), "\n")
}

cat("===========================================================\n")
cat("  Exponential Series System — Simulation Suite\n")
cat("===========================================================\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Experiments:", length(experiments), "\n\n")

t0 <- proc.time()["elapsed"]

for (i in seq_along(experiments)) {
    script <- experiments[i]
    cat("\n===========================================================\n")
    cat("  Running: ", script, "\n")
    cat("===========================================================\n\n")
    source(script, local = TRUE)
}

# --- Export results for vignette ---
cat("\n===========================================================\n")
cat("  Exporting results for vignette\n")
cat("===========================================================\n\n")
source("export_for_vignette.R", local = TRUE)

total <- proc.time()["elapsed"] - t0
cat("\n===========================================================\n")
cat("  All done!\n")
cat("  Total time:", sprintf("%.1f minutes (%.0f seconds)\n", total / 60, total))
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("===========================================================\n")
