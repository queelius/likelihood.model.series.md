# Constructs a likelihood model for `wei_series_homogeneous_md_c1_c2_c3`.

Likelihood model for Weibull series systems with homogeneous shape
parameter and masked component cause of failure with candidate sets that
satisfy conditions C1, C2, and C3.

## Usage

``` r
wei_series_homogeneous_md_c1_c2_c3(
  shape = NULL,
  scales = NULL,
  lifetime = "t",
  indicator = "delta",
  candset = "x"
)
```

## Arguments

- shape:

  common shape parameter for all Weibull component lifetimes

- scales:

  scale parameters for Weibull component lifetimes (optional)

- lifetime:

  column name for system lifetime, defaults to `"t"`

- indicator:

  column name for right-censoring indicator, defaults to `"delta"`.
  TRUE/1 = exact failure time, FALSE/0 = right-censored. For backwards
  compatibility, if this column is not present in the data, censoring is
  inferred from empty candidate sets (all FALSE).

- candset:

  column prefix for candidate set indicators, defaults to `"x"`

## Value

likelihood model object with class
`c("wei_series_homogeneous_md_c1_c2_c3", "series_md", "likelihood_model")`

## Details

This is a reduced model where all components share a common shape
parameter k, while retaining individual scale parameters. The parameter
vector is (k, scale_1, ..., scale_m), giving m+1 parameters instead of
2m.

A key property of this model is that the series system lifetime is
itself Weibull distributed with shape k and scale lambda_s =
(sum(scale_j^-k))^-1/k.

This model satisfies the concept of a `likelihood_model` in the
`likelihood.model` package by providing the following methods:

\(1\) `loglik.wei_series_homogeneous_md_c1_c2_c3` (2)
`score.wei_series_homogeneous_md_c1_c2_c3` (3)
`hess_loglik.wei_series_homogeneous_md_c1_c2_c3`

In this likelihood model, masked component data approximately satisfies:

C1: `Pr{K[i] in C[i]} = 1` C2:
`Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]} = Pr{C[i]=c[i] | K[i]=j', T[i]=t[i]}`
for any `j, j' in c[i]`. C3: masking probabilities are independent of
`theta`

## See also

[`wei_series_md_c1_c2_c3()`](https://queelius.github.io/likelihood.model.series.md/reference/wei_series_md_c1_c2_c3.md)
for the full model with heterogeneous shapes

## Examples

``` r
# Create model and fit to data using generic dispatch
model <- wei_series_homogeneous_md_c1_c2_c3()
# solver <- fit(model)
# theta: (shape, scale1, scale2, ...)
# mle <- solver(data, par = c(1.2, 1000, 900, 850))
```
