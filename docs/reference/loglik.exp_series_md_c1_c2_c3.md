# Log-likelihood method for `exp_series_md_c1_c2_c3` model.

Returns a log-likelihood function for an exponential series system with
respect to rate parameters for masked data with candidate sets that
satisfy conditions C1, C2, and C3.

## Usage

``` r
# S3 method for class 'exp_series_md_c1_c2_c3'
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

- `par`: rate parameters of exponential component lifetime distributions

- `lifetime`: system lifetime column name (default from model)

- `indicator`: right-censoring indicator column name (default from
  model)

- `candset`: prefix of Boolean matrix encoding candidate sets (default
  from model)
