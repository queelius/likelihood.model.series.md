# Hessian of log-likelihood method for `wei_series_homogeneous_md_c1_c2_c3`.

Returns the Hessian (second derivative matrix) of the log-likelihood for
a Weibull series system with homogeneous shape. Computed numerically via
the Jacobian of the score.

## Usage

``` r
# S3 method for class 'wei_series_homogeneous_md_c1_c2_c3'
hess_loglik(model, ...)
```

## Arguments

- model:

  the likelihood model object

- ...:

  additional arguments (passed to returned function)

## Value

hessian function that takes the following arguments:

- `df`: masked data frame

- `par`: parameter vector (shape, scale1, scale2, ...)

- `lifetime`: system lifetime column name (default from model)

- `lifetime_upper`: interval upper bound column name (default from
  model)

- `omega`: observation type column name (default from model)

- `candset`: prefix of Boolean matrix encoding candidate sets
