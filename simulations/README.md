# Simulation Suite: Exponential Series System MLE

Standalone, resumable simulation framework for studying the maximum likelihood
estimator of component failure rates in exponential series systems with masked
data (conditions C1, C2, C3).

## Quick Start

```bash
cd simulations

# Run everything (~30-60 min depending on hardware)
Rscript run_all.R

# Run specific experiments (by number)
Rscript run_all.R 1 3

# Resume after interruption (automatic — just re-run)
Rscript run_all.R

# Force fresh start
Rscript run_all.R --reset

# Run a single experiment
Rscript experiments/01_mc_study.R
```

## Experiments

| # | Script | Description | Key Parameters |
|---|--------|-------------|----------------|
| 1 | `01_mc_study.R` | Main Monte Carlo study: bias, variance, MSE, coverage | n=7500, B=200, p=0.3, ~25% cens |
| 2 | `02_masking_sensitivity.R` | Effect of masking probability on accuracy | n=500, B=100, p∈{0,...,0.5} |
| 3 | `03_censoring_sensitivity.R` | Effect of right-censoring rate on accuracy | n=500, B=100, q∈{0.1,...,0.9} |
| 4 | `04_sample_size.R` | Sample size study with coverage analysis | n∈{100,...,5000}, B=100 |
| 5 | `05_joint_masking_censoring.R` | Joint masking × censoring interaction grid | 3×5 grid, B=100 |

## Directory Structure

```
simulations/
├── R/
│   └── sim_helpers.R          # Shared: data gen, fitting, checkpointing, progress
├── experiments/
│   ├── config.R               # All experimental parameters in one place
│   ├── 01_mc_study.R
│   ├── 02_masking_sensitivity.R
│   ├── 03_censoring_sensitivity.R
│   ├── 04_sample_size.R
│   └── 05_joint_masking_censoring.R
├── results/                   # RDS checkpoints (gitignored)
├── run_all.R                  # Orchestrator
├── export_for_vignette.R      # Transform results → vignette format
└── README.md
```

## Resumability

Each experiment saves its state to `results/<name>.rds` every few replications.
If interrupted (Ctrl-C, crash, system restart), re-running the same script
automatically detects completed work and continues from where it stopped.

The checkpoint uses atomic writes (write to tempfile, then rename) to prevent
corruption from mid-write interruptions.

To force a clean re-run of a specific experiment:

```bash
Rscript experiments/01_mc_study.R --reset
```

## Reproducibility

Each replication uses a deterministic seed derived from the base seed (7231)
and the replication index: `set.seed(SEED + offset + b)`. This means:

- Results are identical regardless of when the simulation was interrupted/resumed
- Individual replications can be reproduced in isolation
- Different experiments use non-overlapping seed ranges to avoid correlation

## Configuration

All parameters are centralized in `experiments/config.R`. To modify an
experiment (e.g., increase B or change the parameter grid), edit config.R
and re-run with `--reset`.

## Vignette Integration

After running experiments 1–3 (minimum), export results for the vignette:

```bash
Rscript export_for_vignette.R
```

This creates `vignettes/precomputed_results.rds` in the format expected by
`vignettes/exponential_series.Rmd`. The vignette loads this file when
`run_long = FALSE` (the default) and renders without re-running simulations.

Experiments 4–5 are included in the export under separate keys
(`sample_size_results`, `joint_results`) for future vignette expansion.

## Dependencies

- `maskedcauses` (this package)
- `md.tools` (for `md_encode_matrix`)
- `dplyr` (for `%>%`)

Install with:

```r
devtools::install_deps(dependencies = TRUE)
```
