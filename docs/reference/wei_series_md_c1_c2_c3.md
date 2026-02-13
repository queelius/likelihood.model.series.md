# Constructs a likelihood model for `wei_series_md_c1_c2_c3`.

Likelihood model for Weibull series systems with masked component cause
of failure with candidate sets that satisfy conditions C1, C2, and C3.

## Usage

``` r
wei_series_md_c1_c2_c3(
  shapes = NULL,
  scales = NULL,
  lifetime = "t",
  lifetime_upper = "t_upper",
  omega = "omega",
  candset = "x"
)
```

## Arguments

- shapes:

  shape parameters for Weibull component lifetimes (optional)

- scales:

  scale parameters for Weibull component lifetimes (optional)

- lifetime:

  column name for system lifetime, defaults to `"t"`

- lifetime_upper:

  column name for interval upper bound, defaults to `"t_upper"`. Only
  used for interval-censored observations.

- omega:

  column name for observation type, defaults to `"omega"`. Must contain
  character values: `"exact"`, `"right"`, `"left"`, or `"interval"`.

- candset:

  column prefix for candidate set indicators, defaults to `"x"`

## Value

likelihood model object with class
`c("wei_series_md_c1_c2_c3", "series_md", "likelihood_model")`

## Details

This model satisfies the concept of a `likelihood_model` in the
`likelihood.model` package by providing the following methods:

\(1\) `loglik.wei_series_md_c1_c2_c3` (2) `score.wei_series_md_c1_c2_c3`
(3) `hess_loglik.wei_series_md_c1_c2_c3`

The Weibull series system has 2m parameters: (shape_1, scale_1, ...,
shape_m, scale_m).

In this likelihood model, masked component data approximately satisfies:

C1: `Pr{K[i] in C[i]} = 1` C2:
`Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]} = Pr{C[i]=c[i] | K[i]=j', T[i]=t[i]}`
for any `j, j' in c[i]`. C3: masking probabilities are independent of
`theta`

## Examples

``` r
# Create model and fit to data using generic dispatch
model <- wei_series_md_c1_c2_c3()
# solver <- fit(model)
# theta: (shape1, scale1, shape2, scale2, ...)
# mle <- solver(data, par = c(1, 1000, 1, 1000, 1, 1000))
```
