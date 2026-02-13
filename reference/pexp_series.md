# Cumulative distribution function for exponential series.

Cumulative distribution function for exponential series.

## Usage

``` r
pexp_series(t, rates, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- t:

  Vector of series system lifetimes.

- rates:

  Vector of rate parameters for exponential component lifetimes.

- lower.tail:

  Logical; if TRUE (default), probabilities are P(X \<= x), otherwise,
  P(X \> x).

- log.p:

  Logical; if TRUE, return the log of the cdf.

## Value

The cumulative probabilities evaluated at the specified lifetimes.
