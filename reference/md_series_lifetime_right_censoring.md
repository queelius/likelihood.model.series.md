# Masked data generation for series system lifetime data

Generates right-censored system failure times and right-censoring
indicators for a series system with the given data frame of component
lifetimes.

## Usage

``` r
md_series_lifetime_right_censoring(
  df,
  tau = Inf,
  comp = "t",
  lifetime = "t",
  right_censoring_indicator = "delta"
)
```

## Arguments

- df:

  a data frame with the indicated component lifetimes

- tau:

  vector of right-censoring times, defaults to `Inf` (no right
  censoring)

- comp:

  component lifetime prefix variable name, defaults to `t`, e.g.,
  `t1, t2, t3`.

- lifetime:

  system lifetime variable name, defaults to `t`

- right_censoring_indicator:

  right-censoring indicator variable, defaults to `delta`

## Value

(masked) data frame with masked data as described in the paper
