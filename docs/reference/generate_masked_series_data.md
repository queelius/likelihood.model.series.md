# Generate masked series system data

Shared data generation logic for all rdata methods. Takes pre-generated
component lifetimes and applies system lifetime calculation,
right-censoring, and candidate set generation.

## Usage

``` r
generate_masked_series_data(
  comp_lifetimes,
  n,
  m,
  tau,
  p,
  default_lifetime,
  default_indicator,
  default_candset
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

  right-censoring time

- p:

  masking probability for non-failed components

- default_lifetime:

  column name for system lifetime

- default_indicator:

  column name for censoring indicator

- default_candset:

  column prefix for candidate sets

## Value

data frame with system lifetime, censoring indicator, and candidate sets
