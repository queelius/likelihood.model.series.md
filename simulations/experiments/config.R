# =============================================================================
# Simulation Configuration
# =============================================================================
# Shared parameters for all experiments in the exponential series system
# masked data MLE simulation study.
#
# Usage: source("experiments/config.R")
# =============================================================================

# --- True model parameters ---
THETA  <- c(1.00, 1.10, 0.95, 1.15, 1.10)
M      <- length(THETA)

# --- Optimizer settings ---
THETA0 <- rep(1, M)         # Initial guess for all components
METHOD <- "Nelder-Mead"     # Optimization method for optim()
ALPHA  <- 0.05              # Significance level for CIs

# --- Reproducibility ---
SEED <- 7231                # Base random seed

# --- Checkpointing ---
SAVE_EVERY <- 5L            # Save checkpoint every N replications

# =============================================================================
# Experiment 1: Main Monte Carlo study
# =============================================================================
MC_N <- 7500L    # Sample size
MC_B <- 200L     # Number of replications
MC_P <- 0.3      # Masking probability
MC_Q <- 0.25     # Survival probability (=> ~25% censoring)

# =============================================================================
# Experiments 2 & 3: Sensitivity analyses
# =============================================================================
SENS_N <- 500L   # Sample size for sensitivity studies
SENS_B <- 100L   # Replications per parameter value

# Experiment 2: Masking sensitivity
MASK_P_VALUES <- seq(0, 0.5, by = 0.1)
MASK_Q        <- 0.25    # Fixed censoring level

# Experiment 3: Censoring sensitivity (9 values for smooth curve)
CENS_Q_VALUES <- seq(0.1, 0.9, by = 0.1)  # Survival probabilities
CENS_P        <- 0.2     # Fixed masking probability

# =============================================================================
# Experiment 4: Sample size study (includes coverage analysis)
# =============================================================================
SIZE_N_VALUES <- c(100L, 250L, 500L, 1000L, 2500L, 5000L)
SIZE_B        <- 100L
SIZE_P        <- 0.3
SIZE_Q        <- 0.25

# =============================================================================
# Experiment 5: Joint masking x censoring interaction
# =============================================================================
JOINT_P_VALUES <- c(0, 0.2, 0.4)
JOINT_Q_VALUES <- c(0.1, 0.3, 0.5, 0.7, 0.9)
JOINT_N        <- 500L
JOINT_B        <- 100L
