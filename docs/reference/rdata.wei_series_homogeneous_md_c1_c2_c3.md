# Random data generation for `wei_series_homogeneous_md_c1_c2_c3` model.

Returns a function that generates random masked data from the
homogeneous Weibull series system DGP at a given parameter value.

## Usage

``` r
# S3 method for class 'wei_series_homogeneous_md_c1_c2_c3'
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
