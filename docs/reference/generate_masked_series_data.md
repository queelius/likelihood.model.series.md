# Generate masked series system data

Shared data generation logic for all rdata methods. Takes pre-generated
component lifetimes and applies an observation mechanism, then generates
candidate sets satisfying conditions C1, C2, C3.

## Usage

``` r
generate_masked_series_data(
  comp_lifetimes,
  n,
  m,
  tau,
  p,
  default_lifetime,
  default_omega,
  default_candset,
  default_lifetime_upper = paste0(default_lifetime, "_upper"),
  observe = NULL
)
```

## Arguments

- comp_lifetimes:

  n x m matrix of component lifetimes

- n:

  number of observations

- m:

  number of components

- tau:

  right-censoring time (used when `observe` is NULL)

- p:

  masking probability for non-failed components

- default_lifetime:

  column name for system lifetime

- default_omega:

  column name for observation type

- default_candset:

  column prefix for candidate sets

- default_lifetime_upper:

  column name for interval upper bound

- observe:

  observation functor created by `observe_*` functions. When NULL, uses
  [`observe_right_censor`](https://queelius.github.io/likelihood.model.series.md/reference/observe_right_censor.md)`(tau)`
  for backwards compatibility.

## Value

data frame with system lifetime, observation type, and candidate sets
