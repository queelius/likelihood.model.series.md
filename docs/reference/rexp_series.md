# Random number generation for exponential series.

Generates random variates from an exponential series distribution.

## Usage

``` r
rexp_series(n, rates, keep_latent = FALSE)
```

## Arguments

- n:

  Integer, number of observations to generate.

- rates:

  Vector of rate parameters for exponential component lifetimes.

- keep_latent:

  Logical; if TRUE, returns a matrix with system lifetimes in the first
  column and individual component lifetimes in subsequent columns. If
  FALSE (default), returns only system lifetimes.

## Value

If `keep_latent = FALSE`, a vector of random variates from the
exponential series distribution. If `keep_latent = TRUE`, a matrix with
system lifetime in the first column and component lifetimes in columns 2
through m+1.
