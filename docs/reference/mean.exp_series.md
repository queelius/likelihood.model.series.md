# Mean function for exponential series.

Computes the expected value of a series system with exponentially
distributed component lifetimes. For a series system with component
rates λ₁,...,λₘ, the system lifetime is exponential with rate Σλⱼ, so
E[T](https://rdrr.io/r/base/logical.html) = 1/Σλⱼ.

## Usage

``` r
# S3 method for class 'exp_series'
mean(x, ...)
```

## Arguments

- x:

  An object of class `exp_series` (a vector of rate parameters).

- ...:

  Additional arguments (ignored, for S3 generic compatibility).

## Value

The mean of the exponential series distribution (1/sum of rates).
