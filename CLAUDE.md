# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

This is an R package for maximum likelihood estimation (MLE) of series system parameters from masked component failure data. The package implements likelihood models for systems where:
- Component lifetimes may be censored (right, left, or interval) via composable observation functors
- Component cause of failure is masked via candidate sets (subsets of components that plausibly contain the failed component)
- Masking conditions C1, C2, C3 are satisfied

## Development Commands

```bash
# Install dependencies (including GitHub remotes)
Rscript -e "devtools::install_deps(dependencies = TRUE)"

# Load package for interactive development
Rscript -e "devtools::load_all()"

# Generate documentation from roxygen2 comments
Rscript -e "devtools::document()"

# Run all tests
Rscript -e "devtools::test()"

# Run a specific test file
Rscript -e "devtools::test(filter = 'exp_series_md')"

# Run tests with coverage (overall)
Rscript -e "covr::package_coverage()"

# Find exact uncovered lines
Rscript -e "covr::zero_coverage(covr::package_coverage())"

# Run R CMD check
Rscript -e "devtools::check()"

# Build and deploy pkgdown site (deploys from docs/ on master)
Rscript -e "pkgdown::build_site()"
```

## Architecture

### Likelihood Models

Three likelihood model classes, all using S3 generic dispatch:

| Model | Parameters | File |
|-------|------------|------|
| `exp_series_md_c1_c2_c3` | m rates: `(λ₁,...,λₘ)` | `R/exp_series_md_c1_c2_c3.R` |
| `wei_series_md_c1_c2_c3` | 2m params: `(k₁,β₁,...,kₘ,βₘ)` | `R/wei_series_md_c1_c2_c3.R` |
| `wei_series_homogeneous_md_c1_c2_c3` | m+1 params: `(k,β₁,...,βₘ)` | `R/wei_series_homogeneous_md_c1_c2_c3.R` |

Each model implements S3 methods: `loglik()`, `score()`, `hess_loglik()`, `assumptions()`, `rdata()`, `fit()`,
`ncomponents()`, `component_hazard()`, `conditional_cause_probability()`, `cause_probability()`

**S3 class hierarchy** (multiple inheritance for dispatch):
```
c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")
c("wei_series_md_c1_c2_c3", "series_md", "likelihood_model")
c("wei_series_homogeneous_md_c1_c2_c3", "series_md", "likelihood_model")
```

**Implementation notes:**
- Exponential: analytical loglik, score, and Hessian for all 4 observation types
- Weibull homogeneous: analytical exact/right, closed-form left/interval via w_c weights; `numDeriv::jacobian` for Hessian
- Weibull heterogeneous: analytical exact/right, `stats::integrate` for left/interval; `numDeriv::jacobian` for Hessian
- Weibull MLE uses L-BFGS-B with positivity lower bounds (Nelder-Mead has poor convergence)

**Usage pattern:**
```r
model <- exp_series_md_c1_c2_c3()
ll_fn <- loglik(model)     # Returns function(df, par, ...)
s_fn <- score(model)       # Returns function(df, par, ...)
H_fn <- hess_loglik(model) # Returns function(df, par, ...)
gen_fn <- rdata(model)     # Returns function(theta, n, tau, p, observe, ...)

# Evaluate
ll_value <- ll_fn(df, par = c(0.5, 0.3, 0.2))

# Generate data with observation functor
df <- gen_fn(theta = c(0.5, 0.3, 0.2), n = 100, p = 0.5,
             observe = observe_periodic(delta = 10, tau = 100))
```

### Observation Functors (`R/observe.R`)

Composable functions that separate the DGP from the observation mechanism.
Signature: `function(t_true) -> list(t, omega, t_upper)`.

| Functor | Description |
|---------|-------------|
| `observe_right_censor(tau)` | Exact if t < tau, right-censored otherwise |
| `observe_left_censor(tau)` | Left-censored if t < tau, right-censored otherwise |
| `observe_periodic(delta, tau)` | Interval-censored at inspection grid; right-censored past tau |
| `observe_mixture(..., weights)` | Randomly selects among schemes per observation |

Pass to `rdata()` via `observe` parameter. Default (NULL) uses `observe_right_censor(tau)`.

### Masking Conditions

- **C1**: Failed component is in candidate set with probability 1
- **C2**: Uniform probability for candidate sets conditioned on component cause
- **C3**: Masking probabilities independent of system parameters

### Data Flow

**Modern path (preferred):** `rdata(model)(theta, n, p, observe = ...)` handles everything internally.

**Manual pipeline (legacy, used in exponential vignette tutorial):**
1. Generate component lifetimes matrix
2. `md_encode_matrix(mat, "t")` → data frame with `t1, t2, ...` columns
3. `md_series_lifetime_right_censoring(df, tau)` → adds `t`, `delta`
4. `md_bernoulli_cand_c1_c2_c3(df, p)` → adds `q1, q2, ...`
5. `md_cand_sampler(df)` → adds `x1, x2, ...`
6. Add `omega` column: `ifelse(delta, "exact", "right")`

### Column Naming Conventions

- `t1, t2, ..., tm`: Component lifetimes (latent)
- `t`: System lifetime (observed)
- `omega`: Observation type — `"exact"`, `"right"`, `"left"`, or `"interval"`
- `t_upper`: Upper bound for interval-censored observations
- `x1, x2, ..., xm`: Boolean candidate set indicators (observed)
- `q1, q2, ..., qm`: Candidate set probabilities (latent, used in manual pipeline)

Model constructors accept column name overrides: `lifetime`, `omega`, `candset`, `lifetime_upper`.

### Matrix Utilities (`R/md_compat.R`)

Internal implementations (formerly from md.tools, now zero external dependencies):
- `md_decode_matrix(df, var)`: Extract prefixed columns as matrix (e.g., `x1, x2` → matrix)
- `md_encode_matrix(mat, var)`: Matrix → data frame with prefixed columns (exported)
- `md_mark_latent(md, vars)`: Set `"latent"` attribute on data frame
- `md_boolean_matrix_to_charsets(df, setvar)`: Boolean matrix → `{1, 3}` display (exported)

### Exponential Series Distribution (`R/exp_series.R`)

Standalone R-style distribution functions (not part of likelihood models):
`dexp_series`, `pexp_series`, `qexp_series`, `rexp_series`, `hazard_exp_series`, `surv.exp_series`

### Series System Generics (`R/series_md_generics.R`)

Package-defined S3 generics dispatched on `series_md` class:
- `ncomponents(model)` — number of components m
- `component_hazard(model, j)` — returns `function(t, par)` for j-th component
- `conditional_cause_probability(model, j)` — returns `function(t, par)` for P(K=j|T=t)
- `cause_probability(model, j)` — returns `function(par)` for marginal P(K=j)

### External Dependencies (GitHub Remotes)

- `algebraic.mle`: MLE objects with `$theta.hat`, `$sigma`, `$info`
- `algebraic.dist`: Distribution utilities
- `likelihood.model`: Generics `loglik`, `score`, `hess_loglik`, `assumptions`, `fim`, `observed_info`, `rdata`
- `generics`: `fit` generic

**Note:** `fim` and `observed_info` are re-exported from `likelihood.model` and use
default implementations (Monte Carlo FIM, negated Hessian for observed info).
No model-specific overrides are needed.

### Utility Functions (`R/utils.R`)

- `cum_haz(haz)`: Creates cumulative hazard from hazard function
- `qcomp(p, haz, surv, theta, t_lower, t_upper)`: Quantile via root-finding
- `rcomp(n, haz, surv, theta)`: Random generation via inverse transform
- `extract_model_data(df, ...)`: Shared validation/extraction for all model methods
- `generate_masked_series_data(...)`: Internal data generation with observe functor support

### Key Mathematical Properties

**Exponential series:** System lifetime T = min(T₁,...,Tₘ) ~ Exp(λ = Σλⱼ)

**Homogeneous Weibull:** When all shapes equal k, system lifetime ~ Weibull(k, βₛ) where βₛ = (Σβⱼ⁻ᵏ)⁻¹/ᵏ. Use `wei_series_system_scale(k, scales)` to compute.

**Log-likelihood (C1, C2, C3) — four observation types:**
- Exact: `log(Σⱼ∈Cᵢ hⱼ(tᵢ)) - H(tᵢ)` where H is system cumulative hazard
- Right: `-H(tᵢ)`
- Left: `log(Σⱼ∈Cᵢ ∫₀ᵗ hⱼ(u)R(u)du) - log(1-R(t))` (or closed form for exp/homogeneous Weibull)
- Interval: `log(Σⱼ∈Cᵢ ∫ₐᵇ hⱼ(u)R(u)du) - log(R(a)-R(b))`

## Testing

Tests use `testthat` edition 3. 798 tests, 99.45% coverage. Test files mirror source files:

| Source | Test |
|--------|------|
| `R/exp_series.R` | `tests/testthat/test-exp_series.R` |
| `R/exp_series_md_c1_c2_c3.R` | `tests/testthat/test-exp_series_md_c1_c2_c3.R` |
| `R/wei_series_md_c1_c2_c3.R` | `tests/testthat/test-wei_series_md_c1_c2_c3.R` |
| `R/wei_series_homogeneous_md_c1_c2_c3.R` | `tests/testthat/test-wei_series_homogeneous_md_c1_c2_c3.R` |
| `R/md_series_*.R` | `tests/testthat/test-md_series.R` |
| `R/utils.R` | `tests/testthat/test-utils.R` |
| `R/observe.R` | `tests/testthat/test-observe.R` |
| `R/series_md_generics.R` | `tests/testthat/test-series_md_generics.R` |
| Censoring integration | `tests/testthat/test-censoring.R` |

## Vignettes

All vignettes use precomputed `.rds` files (set `run_long <- TRUE` to regenerate).

| Vignette | Description |
|----------|-------------|
| `exponential_series.Rmd` | Tutorial: MLE, Monte Carlo simulation, sensitivity analysis |
| `weibull_homogeneous_series.Rmd` | Shared-shape Weibull: DFR/CFR/IFR regimes, closed-form left/interval |
| `weibull_series.Rmd` | Heterogeneous Weibull: mixed hazard shapes, numerical integration |
| `censoring_comparison.Rmd` | Cross-model analysis of censoring type effects on estimation quality |

## Deployment

GitHub Pages serves from `docs/` on master:
```bash
Rscript -e "pkgdown::build_site()"
git add docs/ && git commit && git push
```

## Naming Conventions

- Functions: snake_case (e.g., `md_series_lifetime_right_censoring`)
- Distributions: R conventions (`d*`, `p*`, `q*`, `r*`)
- S3 registration: `#' @method generic class` + `#' @export`
- Observation functors: `observe_*` prefix
- Model constructors accept column name overrides: `lifetime`, `omega`, `candset`, `lifetime_upper`
