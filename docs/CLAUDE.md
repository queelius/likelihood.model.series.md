# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package Overview

This is an R package for maximum likelihood estimation (MLE) of series
system parameters from masked component failure data. The package
implements likelihood models for systems where: - Component lifetimes
may be masked (right-censored) - Component cause of failure is masked
via candidate sets (subsets of components that plausibly contain the
failed component) - Masking conditions C1, C2, C3 are satisfied

## Development Commands

``` bash
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

# Run tests with coverage
Rscript -e "covr::package_coverage()"

# Run R CMD check
Rscript -e "devtools::check()"

# Build pkgdown site
Rscript -e "pkgdown::build_site()"
```

## Architecture

### Likelihood Models

Three likelihood model classes, all using S3 generic dispatch:

| Model                                | Parameters                     | File                                     |
|--------------------------------------|--------------------------------|------------------------------------------|
| `exp_series_md_c1_c2_c3`             | m rates: `(λ₁,...,λₘ)`         | `R/exp_series_md_c1_c2_c3.R`             |
| `wei_series_md_c1_c2_c3`             | 2m params: `(k₁,β₁,...,kₘ,βₘ)` | `R/wei_series_md_c1_c2_c3.R`             |
| `wei_series_homogeneous_md_c1_c2_c3` | m+1 params: `(k,β₁,...,βₘ)`    | `R/wei_series_homogeneous_md_c1_c2_c3.R` |

Each model implements S3 methods:
[`loglik()`](https://queelius.github.io/likelihood.model/reference/loglik.html),
[`score()`](https://queelius.github.io/likelihood.model/reference/score.html),
[`hess_loglik()`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html),
[`assumptions()`](https://queelius.github.io/likelihood.model/reference/assumptions.html),
[`rdata()`](https://queelius.github.io/likelihood.model/reference/rdata.html)

**S3 class hierarchy** (multiple inheritance for dispatch):

    c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")
    c("wei_series_md_c1_c2_c3", "series_md", "likelihood_model")
    c("wei_series_homogeneous_md_c1_c2_c3", "series_md", "likelihood_model")

**Implementation note:** The exponential model computes score and
Hessian analytically. The Weibull models compute score analytically but
use
[`numDeriv::jacobian`](https://rdrr.io/pkg/numDeriv/man/jacobian.html)
on the score for the Hessian.

**Usage pattern:**

``` r
model <- exp_series_md_c1_c2_c3()
ll_fn <- loglik(model)     # Returns function(df, par, ...)
s_fn <- score(model)       # Returns function(df, par, ...)
H_fn <- hess_loglik(model) # Returns function(df, par, ...)
gen_fn <- rdata(model)     # Returns function(theta, n, tau, p, ...)

# Evaluate
ll_value <- ll_fn(df, par = c(0.5, 0.3, 0.2))

# Generate simulated masked data
df <- gen_fn(theta = c(0.5, 0.3, 0.2), n = 100, tau = 10, p = 0.5)
```

### Masking Conditions

- **C1**: Failed component is in candidate set with probability 1
- **C2**: Uniform probability for candidate sets conditioned on
  component cause
- **C3**: Masking probabilities independent of system parameters

### Data Flow

1.  Generate component lifetimes (`t1, t2, ...`)
2.  Apply
    [`md_series_lifetime_right_censoring()`](https://queelius.github.io/likelihood.model.series.md/reference/md_series_lifetime_right_censoring.md)
    → adds `t`, `delta`
3.  Apply
    [`md_bernoulli_cand_c1_c2_c3()`](https://queelius.github.io/likelihood.model.series.md/reference/md_bernoulli_cand_c1_c2_c3.md)
    → adds `q1, q2, ...`
4.  Apply
    [`md_cand_sampler()`](https://queelius.github.io/likelihood.model.series.md/reference/md_cand_sampler.md)
    → adds `x1, x2, ...`
5.  Fit model using `t` and `x1, x2, ...` columns

### Column Naming Conventions

- `t1, t2, ..., tm`: Component lifetimes (latent)
- `t`: System lifetime (observed)
- `delta`: Right-censoring indicator (TRUE = exact, FALSE = censored)
- `q1, q2, ..., qm`: Candidate set probabilities (latent)
- `x1, x2, ..., xm`: Boolean candidate set indicators (observed)

### Exponential Series Distribution (`R/exp_series.R`)

Standalone R-style distribution functions (not part of likelihood
models): `dexp_series`, `pexp_series`, `qexp_series`, `rexp_series`,
`hazard_exp_series`, `surv.exp_series`

### External Dependencies (GitHub Remotes)

- `md.tools`:
  [`md_decode_matrix()`](https://queelius.github.io/likelihood.model.series.md/reference/md_decode_matrix.md),
  [`md_encode_matrix()`](https://queelius.github.io/likelihood.model.series.md/reference/md_encode_matrix.md),
  [`md_mark_latent()`](https://queelius.github.io/likelihood.model.series.md/reference/md_mark_latent.md)
- `algebraic.mle`: MLE objects with `$theta.hat`, `$sigma`, `$info`
- `algebraic.dist`: Distribution utilities
- `likelihood.model`: Generics `loglik`, `score`, `hess_loglik`,
  `assumptions`, `fim`, `observed_info`, `rdata`
- `generics`: `fit` generic

**Note:** `fim` and `observed_info` are re-exported from
`likelihood.model` and use default implementations (Monte Carlo FIM,
negated Hessian for observed info). No model-specific overrides are
needed.

### Key Mathematical Properties

**Exponential series:** System lifetime T = min(T₁,…,Tₘ) ~ Exp(λ = Σλⱼ)

**Homogeneous Weibull:** When all shapes equal k, system lifetime ~
Weibull(k, βₛ) where βₛ = (Σβⱼ⁻ᵏ)⁻¹/ᵏ. Use
`wei_series_system_scale(k, scales)` to compute.

**Log-likelihood (C1, C2, C3):**

    ℓ(θ) = -Σᵢ [survival term] + Σᵢ:δᵢ=1 log(Σⱼ∈Cᵢ hⱼ(tᵢ))

### Utility Functions (`R/utils.R`)

- `cum_haz(haz)`: Creates cumulative hazard from hazard function
- `qcomp(p, haz, surv, theta, t_lower, t_upper)`: Quantile via
  root-finding
- `rcomp(n, haz, surv, theta)`: Random generation via inverse transform

## Testing

Tests use `testthat` edition 3. Test files mirror source files:

| Source                                   | Test                                                       |
|------------------------------------------|------------------------------------------------------------|
| `R/exp_series.R`                         | `tests/testthat/test-exp_series.R`                         |
| `R/exp_series_md_c1_c2_c3.R`             | `tests/testthat/test-exp_series_md_c1_c2_c3.R`             |
| `R/wei_series_md_c1_c2_c3.R`             | `tests/testthat/test-wei_series_md_c1_c2_c3.R`             |
| `R/wei_series_homogeneous_md_c1_c2_c3.R` | `tests/testthat/test-wei_series_homogeneous_md_c1_c2_c3.R` |
| `R/md_series_*.R`                        | `tests/testthat/test-md_series.R`                          |
| `R/utils.R`                              | `tests/testthat/test-utils.R`                              |

## Vignettes

| Vignette                           | Description                                                 |
|------------------------------------|-------------------------------------------------------------|
| `vignettes/exponential_series.Rmd` | Tutorial: MLE, Monte Carlo simulation, sensitivity analysis |

## Naming Conventions

- Functions: snake_case (e.g., `md_series_lifetime_right_censoring`)
- Distributions: R conventions (`d*`, `p*`, `q*`, `r*`)
- S3 registration: `#' @method generic class` + `#' @export`
- Model constructors accept column name overrides: `lifetime`,
  `indicator`, `candset` params

## API Notes

**Unified Censoring API:** All models accept an `indicator` parameter
(default: `"delta"`) specifying the column for right-censoring status.
For backwards compatibility, if this column is absent, censoring is
inferred from empty candidate sets (all FALSE).
