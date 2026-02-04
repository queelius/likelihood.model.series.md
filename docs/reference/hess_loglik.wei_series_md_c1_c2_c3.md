# Hessian of log-likelihood method for `wei_series_md_c1_c2_c3` model.

Returns the Hessian (second derivative matrix) of the log-likelihood for
a Weibull series system. Computed numerically via the Jacobian of the
score.

## Usage

``` r
# S3 method for class 'wei_series_md_c1_c2_c3'
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

- `par`: parameter vector (shape1, scale1, shape2, scale2, ...)

- `lifetime`: system lifetime column name (default from model)

- `indicator`: right-censoring indicator column name (default from
  model)

- `candset`: prefix of Boolean matrix encoding candidate sets (default
  from model)
