# Quantile function for exponential series.

Quantile function for exponential series.

## Usage

``` r
qexp_series(p, rates, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- p:

  Vector of quantiles.

- rates:

  Vector of rate parameters for exponential component lifetimes.

- lower.tail:

  Logical, if TRUE (default), probabilities are P(X \<= x), otherwise,
  P(X \> x).

- log.p:

  Logical, if TRUE, vector of probabilities `p` are given as `log(p)`.

## Value

Quantiles corresponding to the given probabilities `p`.
