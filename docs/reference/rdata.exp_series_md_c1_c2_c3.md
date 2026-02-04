# Random data generation for `exp_series_md_c1_c2_c3` model.

Returns a function that generates random masked data from the
exponential series system DGP at a given parameter value.

## Usage

``` r
# S3 method for class 'exp_series_md_c1_c2_c3'
rdata(model, ...)
```

## Arguments

- model:

  the likelihood model object

- ...:

  additional arguments (passed to returned function)

## Value

function that takes (theta, n, tau, p, ...) and returns a data frame
with columns: t, delta, x1, x2, ..., xm
