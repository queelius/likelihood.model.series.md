# Score method for `wei_series_md_c1_c2_c3` model.

Returns a score (gradient) function for a Weibull series system with
respect to parameter vector (shape_1, scale_1, ..., shape_m, scale_m)
for masked data with candidate sets that satisfy conditions C1, C2, and
C3.

## Usage

``` r
# S3 method for class 'wei_series_md_c1_c2_c3'
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

- `par`: parameter vector (shape1, scale1, shape2, scale2, ...)

- `lifetime`: system lifetime column name (default from model)

- `lifetime_upper`: interval upper bound column name (default from
  model)

- `omega`: observation type column name (default from model)

- `candset`: prefix of Boolean matrix encoding candidate sets

## Details

Uses a hybrid approach: analytical score for exact and right-censored
observations, numerical gradient (via numDeriv) for left-censored and
interval-censored observations.
