# Survival function for exponential series.

Survival function for exponential series.

## Usage

``` r
surv.exp_series(t, rates, log.p = FALSE)
```

## Arguments

- t:

  Vector of series system lifetimes.

- rates:

  Vector of rate parameters for exponential component lifetimes.

- log.p:

  Logical; if TRUE, return the log of the survival function.

## Value

The survival function evaluated at the specified lifetimes.
