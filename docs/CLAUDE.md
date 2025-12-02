# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package Overview

This is an R package for maximum likelihood estimation (MLE) of series
system parameters from masked component failure data. The package
implements likelihood models for systems where: - Component lifetimes
may be masked (right-censored, left-censored, interval-censored) -
Component cause of failure is masked via candidate sets (subsets of
components that plausibly contain the failed component) - Various
masking conditions (C1, C2, C3) can be satisfied or relaxed

## Development Commands

### Building and Installing

``` bash
# Install package dependencies (including GitHub remotes)
Rscript -e "devtools::install_deps(dependencies = TRUE)"

# Build the package
R CMD build .

# Install the package locally
Rscript -e "devtools::install()"

# Load package for interactive development
Rscript -e "devtools::load_all()"
```

### Testing

``` bash
# Run all tests
Rscript -e "devtools::test()"

# Run tests with coverage
Rscript -e "covr::package_coverage()"
```

**Note:** Currently, this package does not have a `tests/` directory.
When adding tests, create the structure:

``` bash
mkdir -p tests/testthat
echo 'library(testthat)\nlibrary(likelihood.model.series.md)\ntest_check("likelihood.model.series.md")' > tests/testthat.R
```

### Documentation

``` bash
# Generate documentation from roxygen2 comments
Rscript -e "devtools::document()"

# Build pkgdown site
Rscript -e "pkgdown::build_site()"

# Preview pkgdown site locally
Rscript -e "pkgdown::preview_site()"
```

### Checking Package Health

``` bash
# Run R CMD check
Rscript -e "devtools::check()"

# Check examples
Rscript -e "devtools::run_examples()"
```

## Architecture

### Core Likelihood Model Framework

The package is built on the `likelihood.model` framework, which
provides: - Generic likelihood contribution models via
`likelihood_contr_model` - Methods for MLE fitting, bootstrapping,
confidence intervals - Integration with `algebraic.mle` for MLE objects
with methods like [`predict()`](https://rdrr.io/r/stats/predict.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html)

### Masking Conditions

Models are characterized by three conditions: - **C1**: Failed component
is in candidate set with probability 1 - **C2**: Uniform probability for
candidate sets conditioned on component cause - **C3**: Masking
probabilities independent of system parameters

### Key Components

#### 1. Series System Distributions (`R/exp_series.R`)

Exponential series system distribution functions: -
[`rexp_series()`](https://queelius.github.io/likelihood.md.series.systems/reference/rexp_series.md):
Random generation (with optional `keep_latent` parameter to retain
component lifetimes) -
[`dexp_series()`](https://queelius.github.io/likelihood.md.series.systems/reference/dexp_series.md),
[`pexp_series()`](https://queelius.github.io/likelihood.md.series.systems/reference/pexp_series.md),
[`qexp_series()`](https://queelius.github.io/likelihood.md.series.systems/reference/qexp_series.md):
PDF, CDF, quantile functions -
[`hazard_exp_series()`](https://queelius.github.io/likelihood.md.series.systems/reference/hazard_exp_series.md),
[`surv.exp_series()`](https://queelius.github.io/likelihood.md.series.systems/reference/surv.exp_series.md):
Hazard and survival functions -
[`mean.exp_series()`](https://queelius.github.io/likelihood.md.series.systems/reference/mean.exp_series.md):
Mean function for exponential series

**Key implementation detail:** For a series system with exponential
components with rates `λ₁,...,λₘ`, the system lifetime is exponentially
distributed with rate `λ = Σλⱼ`.

#### 2. Likelihood Models (`R/exp_series_md_c1_c2_c3.R`)

Implements likelihood model for exponential series with masked data
satisfying C1, C2, C3: -
[`exp_series_md_c1_c2_c3()`](https://queelius.github.io/likelihood.md.series.systems/reference/exp_series_md_c1_c2_c3.md):
Constructor for likelihood model object (returns S3 object with class
vector `c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")`) -
S3 methods:
[`loglik()`](https://queelius.github.io/likelihood.model/reference/loglik.html),
[`score()`](https://queelius.github.io/likelihood.model/reference/score.html),
[`hess_loglik()`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html),
[`assumptions()`](https://queelius.github.io/likelihood.model/reference/assumptions.html)
dispatch to the model

These functions work with data frames containing: - `t` (or custom):
System lifetime column - `x1, x2, ..., xm` (or custom prefix): Boolean
candidate set indicators

**Generic Dispatch Pattern:** The package uses S3 generic dispatch.
Generics from `likelihood.model` (`loglik`, `score`, `hess_loglik`,
`assumptions`, `fim`) and `generics` (`fit`) are re-exported:

``` r
model <- exp_series_md_c1_c2_c3()

# Get function generators via generic dispatch
ll_fn <- loglik(model)   # Returns function(df, par, ...)
s_fn <- score(model)     # Returns function(df, par, ...)
H_fn <- hess_loglik(model)  # Returns function(df, par, ...)

# Evaluate directly
ll_value <- ll_fn(df, par = c(0.5, 0.3, 0.2))

# Or fit using the generic fit() function
solver <- fit(model)
mle_result <- solver(df, par = c(1, 1, 1))
cat("Estimates:", mle_result$theta.hat)
```

#### 3. Masked Data Generation (`R/md_series_lifetime.R`, `R/md_series_masked_component_cause.R`)

Functions for generating and working with masked data:

**Lifetime Masking:** -
[`md_series_lifetime_right_censoring()`](https://queelius.github.io/likelihood.md.series.systems/reference/md_series_lifetime_right_censoring.md):
Applies right-censoring to series system lifetimes - Reads component
lifetimes from columns `t1, t2, ...` (or custom prefix) - Adds `t`
(system lifetime = min of component lifetimes or censoring time `tau`) -
Adds `delta` (right-censoring indicator: TRUE = exact, FALSE =
right-censored) - Marks component lifetimes as latent using
`md_mark_latent()`

**Component Cause Masking:** -
[`md_bernoulli_cand_c1_c2_c3()`](https://queelius.github.io/likelihood.md.series.systems/reference/md_bernoulli_cand_c1_c2_c3.md):
Generates Bernoulli candidate set probabilities satisfying C1, C2, C3 -
Reads component lifetimes from `t1, t2, ...` columns - Creates
probability columns `q1, q2, ...` where `qⱼ = 1` for failed component,
`qⱼ = p` otherwise - If system is right-censored (`delta = FALSE`), sets
all probabilities to 0 - Returns data frame with added `q` columns
(marked as latent)

- [`md_cand_sampler()`](https://queelius.github.io/likelihood.md.series.systems/reference/md_cand_sampler.md):
  Samples candidate sets from component probabilities
  - Reads probability columns `q1, q2, ...`
  - Samples Boolean candidate set indicators `x1, x2, ...` where
    `xⱼ ~ Bernoulli(qⱼ)`
  - Returns data frame with added `x` columns

### Data Flow Pattern

1.  Generate/load component lifetime data (`t1, t2, ...`)
2.  Apply system lifetime masking → add `t`, `delta` columns via
    [`md_series_lifetime_right_censoring()`](https://queelius.github.io/likelihood.md.series.systems/reference/md_series_lifetime_right_censoring.md)
3.  Generate candidate set probabilities → add `q1, q2, ...` columns via
    [`md_bernoulli_cand_c1_c2_c3()`](https://queelius.github.io/likelihood.md.series.systems/reference/md_bernoulli_cand_c1_c2_c3.md)
4.  Sample candidate sets → add `x1, x2, ...` columns via
    [`md_cand_sampler()`](https://queelius.github.io/likelihood.md.series.systems/reference/md_cand_sampler.md)
5.  Fit likelihood model using `x` columns and `t` column

**Example:**

``` r
library(likelihood.model.series.md)

# Step 1: Generate component lifetimes
df <- data.frame(
  t1 = rexp(100, 0.5),
  t2 = rexp(100, 0.3),
  t3 = rexp(100, 0.2)
)

# Step 2: Apply right-censoring
df <- md_series_lifetime_right_censoring(df, tau = 10)

# Step 3: Generate candidate set probabilities
df <- md_bernoulli_cand_c1_c2_c3(df, p = 0.3)

# Step 4: Sample candidate sets
df <- md_cand_sampler(df)

# Step 5: Fit likelihood model using generic dispatch
model <- exp_series_md_c1_c2_c3()
solver <- fit(model)
mle <- solver(df, par = c(1, 1, 1))

# Access results
mle$theta.hat      # MLE estimates
mle$sigma          # Standard errors
mle$info           # Fisher information matrix
print(mle)         # Pretty-printed summary with CIs
```

### Column Naming Conventions

- `t1, t2, ..., tm`: Component lifetimes (often latent/unobserved)
- `t`: Series system lifetime (minimum of component lifetimes or
  censored time)
- `delta`: Right-censoring indicator (TRUE = exact, FALSE =
  right-censored)
- `q1, q2, ..., qm`: Component inclusion probabilities for candidate
  sets (latent)
- `x1, x2, ..., xm`: Boolean candidate set indicators (observed)

**Matrix encoding/decoding:** Use
[`md.tools::md_decode_matrix()`](https://queelius.github.io/md.tools/reference/md_decode_matrix.html)
and
[`md.tools::md_encode_matrix()`](https://queelius.github.io/md.tools/reference/md_encode_matrix.html)
to convert between wide format (columns `prefix1, prefix2, ...`) and
matrix format.

### External Dependencies

The package depends on several related packages by the same author
(installed via GitHub remotes in DESCRIPTION): - `md.tools`: Utilities
for encoding/decoding masked data matrices (`md_decode_matrix`,
`md_encode_matrix`, `md_mark_latent`) - `algebraic.mle`: MLE objects
with rich method support (access estimates via `$theta.hat`, `$sigma`,
`$info`) - `algebraic.dist`: Distribution algebra utilities -
`likelihood.model`: General likelihood model framework providing generic
functions (`loglik`, `score`, `hess_loglik`, `assumptions`, `fim`) and
`fit.likelihood_model` - `generics`: Provides the `fit` generic
(re-exported by this package)

**Re-exported Generics:** This package re-exports `loglik`, `score`,
`hess_loglik`, `assumptions`, `fim` from `likelihood.model` and `fit`
from `generics`, so users only need to load this package.

### Utility Functions (`R/utils.R`)

- [`cum_haz()`](https://queelius.github.io/likelihood.md.series.systems/reference/cum_haz.md):
  Creates cumulative hazard function from hazard function
- `qcomp()`: Quantile function for component with custom hazard/survival
- `rcomp()`: Random generation for component with custom hazard/survival

## Naming Conventions

- Functions use snake_case (e.g., `md_series_lifetime_right_censoring`)
- Distribution functions follow R conventions: `d*` (density), `p*`
  (CDF), `q*` (quantile), `r*` (random)
- Model constructors return S3 objects with class vectors like
  `c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")`
- S3 methods are registered properly using `#' @method generic class`
  and `#' @export` roxygen2 tags
- Column prefixes: `t` (time), `q` (probability), `x` (candidate set
  boolean), `delta` (censoring indicator)

## Documentation Standards

- All exported functions use roxygen2 documentation
- Main vignette: `vignettes/exponential_series.Rmd` (comprehensive
  mathematical treatment)
- Package website:
  <https://queelius.github.io/likelihood.model.series.md/>
- Documentation is auto-generated via pkgdown (GitHub Actions workflow
  on push to master/main)
- Use `@importFrom` for specific imports rather than importing entire
  packages

## Development Notes

### Mathematical Foundation

The exponential series system has the property that if components have
rates `λ₁,...,λₘ`, the system lifetime `T = min(T₁,...,Tₘ)` is
exponentially distributed with rate `λ = Σλⱼ`. This simplifies many
calculations.

### Log-Likelihood Formula (C1, C2, C3)

For masked data satisfying C1, C2, C3:

    ℓ(θ|data) = -Σᵢ tᵢ · Σⱼ θⱼ + Σᵢ log(Σⱼ∈Cᵢ θⱼ)

where `Cᵢ` is the candidate set for observation `i`, and `tᵢ` is the
system lifetime.

### Score and Hessian

The score (gradient) and Hessian are computed analytically in
`R/exp_series_md_c1_c2_c3.R`. This allows for efficient Newton-Raphson
optimization.

### Future Extensions

The `.fixing/` directory contains experimental code for: - Weibull
series systems (`wei_series.R`) - Alternative candidate set models
(`cand-set-models.R`) - Sensitivity analysis (`sen-ana.Rmd`)

These are not currently exported or documented.
