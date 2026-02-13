# likelihood.model.series.md

Maximum likelihood estimation for series system parameters from masked
component failure data. Supports exact, right-censored, left-censored,
and interval-censored observations with composable observation schemes.

## Installation

From [r-universe](https://queelius.r-universe.dev):

``` r
install.packages("likelihood.model.series.md",
                 repos = "https://queelius.r-universe.dev")
```

Or from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("queelius/likelihood.model.series.md")
```

## Models

Three likelihood models, all using S3 generic dispatch:

| Model                                | Parameters              | Description                                                                                    |
|--------------------------------------|-------------------------|------------------------------------------------------------------------------------------------|
| `exp_series_md_c1_c2_c3`             | m rates: (λ₁, …, λₘ)    | Exponential components. Closed-form loglik, score, and Hessian for all four observation types. |
| `wei_series_homogeneous_md_c1_c2_c3` | m+1: (k, β₁, …, βₘ)     | Weibull with shared shape. Closed-form loglik; hybrid analytical/numerical score.              |
| `wei_series_md_c1_c2_c3`             | 2m: (k₁, β₁, …, kₘ, βₘ) | Weibull with per-component shapes. Numerical integration for left/interval types.              |

Each model implements S3 methods for:
[`loglik()`](https://queelius.github.io/likelihood.model/reference/loglik.html),
[`score()`](https://queelius.github.io/likelihood.model/reference/score.html),
[`hess_loglik()`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html),
[`rdata()`](https://queelius.github.io/likelihood.model/reference/rdata.html),
[`assumptions()`](https://queelius.github.io/likelihood.model/reference/assumptions.html),
[`component_hazard()`](https://queelius.github.io/likelihood.model.series.md/reference/component_hazard.md),
[`cause_probability()`](https://queelius.github.io/likelihood.model.series.md/reference/cause_probability.md),
[`conditional_cause_probability()`](https://queelius.github.io/likelihood.model.series.md/reference/conditional_cause_probability.md),
[`ncomponents()`](https://queelius.github.io/likelihood.model.series.md/reference/ncomponents.md)

## Quick Start

``` r
library(likelihood.model.series.md)

# Create model and generate data
model <- exp_series_md_c1_c2_c3()
gen <- rdata(model)
df <- gen(theta = c(1.0, 1.1, 0.95), n = 500, tau = 3, p = 0.3)

# Fit via MLE
solver <- fit(model)
mle <- solver(df, par = c(1, 1, 1))
mle$theta.hat   # estimated rates
mle$sigma       # asymptotic standard errors

# Evaluate likelihood functions at the MLE
loglik(model)(df, par = mle$theta.hat)
score(model)(df, par = mle$theta.hat)
hess_loglik(model)(df, par = mle$theta.hat)
```

## Observation Types

The package supports four observation types via composable **observation
functors** that separate the observation mechanism from the
data-generating process:

| Type              | Column `omega` | Description                               |
|-------------------|----------------|-------------------------------------------|
| Exact             | `"exact"`      | Failure time observed precisely           |
| Right-censored    | `"right"`      | System survived past observation time     |
| Left-censored     | `"left"`       | System failed before inspection time      |
| Interval-censored | `"interval"`   | Failure occurred in window `[t, t_upper]` |

``` r
# Right-censoring at tau (default, backwards-compatible)
df <- gen(theta, n = 500, tau = 5, p = 0.3)

# Periodic inspection every 0.5 time units
df <- gen(theta, n = 500, p = 0.3,
          observe = observe_periodic(delta = 0.5, tau = 5))

# Mixed monitoring: 60% continuous, 20% single inspection, 20% periodic
df <- gen(theta, n = 500, p = 0.3,
          observe = observe_mixture(
            observe_right_censor(tau = 5),
            observe_left_censor(tau = 3),
            observe_periodic(delta = 0.5, tau = 5),
            weights = c(0.6, 0.2, 0.2)
          ))
```

## Masking Conditions

The candidate set models satisfy three conditions (C1, C2, C3) that
yield a reduced likelihood depending only on observed data:

- **C1**: The failed component is always in the candidate set
- **C2**: Candidate set probabilities are symmetric across components in
  the set
- **C3**: Masking probabilities are independent of system parameters

## Data Format

| Column            | Description                                                       |
|-------------------|-------------------------------------------------------------------|
| `t`               | System lifetime (observed)                                        |
| `omega`           | Observation type: `"exact"`, `"right"`, `"left"`, or `"interval"` |
| `t_upper`         | Upper bound for interval-censored observations (NA otherwise)     |
| `x1, x2, ..., xm` | Boolean candidate set indicators                                  |

For backwards compatibility, a `delta` column (TRUE/FALSE) is also
accepted in place of `omega` for exact/right-censored data.

## Vignettes

Four tutorial vignettes with worked examples and Monte Carlo studies:

``` r
vignette("exponential_series", package = "likelihood.model.series.md")
vignette("weibull_homogeneous_series", package = "likelihood.model.series.md")
vignette("weibull_series", package = "likelihood.model.series.md")
vignette("censoring_comparison", package = "likelihood.model.series.md")
```

Browse online at
<https://queelius.github.io/likelihood.model.series.md/articles/>.

## Related Packages

- [likelihood.model](https://github.com/queelius/likelihood.model) –
  generics for likelihood models (`loglik`, `score`, `fit`, etc.)
- [algebraic.mle](https://github.com/queelius/algebraic.mle) – MLE
  objects with inference methods
- [algebraic.dist](https://github.com/queelius/algebraic.dist) –
  distribution utilities
- [md.tools](https://github.com/queelius/md.tools) – masked data
  encoding/decoding

## License

GPL (\>= 3)
