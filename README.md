# likelihood.model.series.md

<!-- badges: start -->
[![R-universe](https://queelius.r-universe.dev/badges/likelihood.model.series.md)](https://queelius.r-universe.dev/likelihood.model.series.md)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Maximum likelihood estimation for series system parameters from masked
component failure data. Implements likelihood models for systems where
component lifetimes may be right-censored and the component cause of
failure is identified only up to a candidate set.

## Installation

From [r-universe](https://queelius.r-universe.dev):

```r
install.packages("likelihood.model.series.md",
                 repos = "https://queelius.r-universe.dev")
```

Or from GitHub:

```r
# install.packages("devtools")
devtools::install_github("queelius/likelihood.model.series.md")
```

## Models

Three likelihood models are provided, all using S3 generic dispatch:

| Model | Parameters | Description |
|-------|-----------|-------------|
| `exp_series_md_c1_c2_c3` | m rates (λ₁, ..., λₘ) | Exponential component lifetimes |
| `wei_series_md_c1_c2_c3` | 2m params (k₁, β₁, ..., kₘ, βₘ) | Weibull with per-component shapes |
| `wei_series_homogeneous_md_c1_c2_c3` | m+1 params (k, β₁, ..., βₘ) | Weibull with shared shape |

Each model implements: `loglik()`, `score()`, `hess_loglik()`, `rdata()`, `assumptions()`

## Masking Conditions

The candidate set models satisfy three conditions that yield a reduced
likelihood depending only on observed data:

- **C1**: The failed component is in the candidate set with probability 1
- **C2**: Given that the failed component is in the candidate set, the
  probability of observing any particular candidate set is the same
  regardless of which component in the set actually failed
- **C3**: Masking probabilities are independent of the system parameters θ

## Usage

```r
library(likelihood.model.series.md)

# --- Exponential model ---
model <- exp_series_md_c1_c2_c3()

# Generate masked series data
gen <- rdata(model)
df <- gen(theta = c(1.0, 1.1, 0.95), n = 500, tau = 3, p = 0.3)

# Fit via MLE
solver <- fit(model)
mle <- solver(df, par = c(1, 1, 1))
mle$theta.hat   # estimated rates
mle$sigma       # asymptotic standard errors

# Log-likelihood, score, and Hessian
ll <- loglik(model)
s  <- score(model)
H  <- hess_loglik(model)

ll(df, par = mle$theta.hat)
s(df, par = mle$theta.hat)
H(df, par = mle$theta.hat)

# --- Weibull model (shared shape) ---
wmodel <- wei_series_homogeneous_md_c1_c2_c3()
wsolver <- fit(wmodel)
# theta: (shape, scale1, scale2, ...)
# wmle <- wsolver(data, par = c(1.2, 1000, 900, 850))
```

## Data Format

The data frame passed to model functions uses these column conventions:

| Column | Description |
|--------|-------------|
| `t` | System lifetime (observed) |
| `delta` | Right-censoring indicator (TRUE = exact, FALSE = censored) |
| `x1, x2, ..., xm` | Boolean candidate set indicators |

Column names are configurable via the `lifetime`, `indicator`, and `candset`
parameters when constructing a model.

## Simulation Framework

The `simulations/` directory contains a standalone, resumable Monte Carlo
framework for studying MLE properties under varying masking, censoring, and
sample size conditions. Five experiments are included:

1. Main MC study (bias, variance, MSE, coverage)
2. Masking probability sensitivity
3. Censoring rate sensitivity
4. Sample size scaling
5. Joint masking × censoring interaction

See [`simulations/README.md`](simulations/README.md) for details. Results
can be exported for the package vignette via `simulations/export_for_vignette.R`.

## Vignette

The package includes a tutorial vignette with worked examples and
data-driven analysis of MLE properties:

```r
vignette("exponential_series", package = "likelihood.model.series.md")
```

## Related Packages

- [likelihood.model](https://github.com/queelius/likelihood.model) — generics for likelihood models (`loglik`, `score`, `fit`, etc.)
- [algebraic.mle](https://github.com/queelius/algebraic.mle) — MLE objects with inference methods
- [algebraic.dist](https://github.com/queelius/algebraic.dist) — distribution utilities
- [md.tools](https://github.com/queelius/md.tools) — masked data encoding/decoding

## License

GPL (>= 3)
