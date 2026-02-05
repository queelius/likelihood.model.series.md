# Score method for `wei_series_homogeneous_md_c1_c2_c3` model.

Returns a score (gradient) function for a Weibull series system with
homogeneous shape parameter. The parameter vector is (k, scale_1, ...,
scale_m) for masked data with candidate sets that satisfy conditions C1,
C2, and C3.

## Usage

``` r
# S3 method for class 'wei_series_homogeneous_md_c1_c2_c3'
score(model, ...)
```

## Arguments

- model:

  the likelihood model object

- ...:

  additional arguments (passed to returned function)

## Value

score function that takes the following arguments:

- `df`: masked data frame

- `par`: parameter vector (shape, scale1, scale2, ...)

- `lifetime`: system lifetime column name (default from model)

- `indicator`: right-censoring indicator column name (default from
  model)

- `candset`: prefix of Boolean matrix encoding candidate sets
