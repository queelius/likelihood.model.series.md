# Score method for `exp_series_md_c1_c2_c3` model.

Returns a score (gradient) function for an exponential series system
with respect to parameter `theta` for masked component failure with
candidate sets that satisfy conditions C1, C2, and C3.

## Usage

``` r
# S3 method for class 'exp_series_md_c1_c2_c3'
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

- `par`: rate parameters of exponential component lifetime distributions

- `lifetime`: system lifetime column name (default from model)

- `lifetime_upper`: interval upper bound column name (default from
  model)

- `omega`: observation type column name (default from model)

- `candset`: prefix of Boolean matrix encoding candidate sets

## Details

All four observation types (exact, right, left, interval) have
closed-form score contributions.
