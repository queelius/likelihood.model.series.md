# Hazard function for exponential series.

Hazard function for exponential series.

## Usage

``` r
hazard_exp_series(t, rates, log.p = FALSE)
```

## Arguments

- t:

  Vector of series system lifetimes.

- rates:

  Vector of rate parameters for exponential component lifetimes.

- log.p:

  Logical; if TRUE, return the log of the hazard function.

## Value

The hazard function evaluated at the specified lifetimes.
