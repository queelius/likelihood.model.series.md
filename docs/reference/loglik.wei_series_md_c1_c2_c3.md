# Log-likelihood method for `wei_series_md_c1_c2_c3` model.

Returns a log-likelihood function for a Weibull series system with
respect to parameter vector (shape_1, scale_1, ..., shape_m, scale_m)
for masked data with candidate sets that satisfy conditions C1, C2, and
C3.

## Usage

``` r
# S3 method for class 'wei_series_md_c1_c2_c3'
loglik(model, ...)
```

## Arguments

- model:

  the likelihood model object

- ...:

  additional arguments (passed to returned function)

## Value

log-likelihood function that takes the following arguments:

- `df`: masked data frame

- `par`: parameter vector (shape1, scale1, shape2, scale2, ...)

- `lifetime`: system lifetime column name (default from model)

- `indicator`: right-censoring indicator column name (default from
  model)

- `candset`: prefix of Boolean matrix encoding candidate sets
