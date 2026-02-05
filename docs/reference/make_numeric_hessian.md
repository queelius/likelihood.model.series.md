# Create numeric Hessian function from score function via Jacobian

Factory function for creating Hessian functions using numerical
differentiation of the score function. Used by Weibull models which
compute score analytically but need numerical Hessian.

## Usage

``` r
make_numeric_hessian(score_fn, model)
```

## Arguments

- score_fn:

  score function returned by score.\* method

- model:

  likelihood model object (for extracting defaults)

## Value

function(df, par, ...) that computes the Hessian matrix
