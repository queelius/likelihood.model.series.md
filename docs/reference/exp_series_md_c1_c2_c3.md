# Constructs a likelihood model for `exp_series_md_c1_c2_c3`.

Likelihood model for exponential series systems with masked component
cause of failure with candidate sets that satisfy conditions C1, C2, and
C3, described below.

## Usage

``` r
exp_series_md_c1_c2_c3(
  rates = NULL,
  lifetime = "t",
  lifetime_upper = "t_upper",
  omega = "omega",
  candset = "x"
)
```

## Arguments

- rates:

  rate parameters for exponential component lifetimes (optional, used as
  initial values for MLE if provided)

- lifetime:

  column name for system lifetime, defaults to `"t"`

- lifetime_upper:

  column name for interval upper bound, defaults to `"t_upper"`. Only
  used for interval-censored observations.

- omega:

  column name for observation type, defaults to `"omega"`. Must contain
  character values: `"exact"` (failure at t), `"right"` (right-censored
  at t), `"left"` (left-censored: failed before t), or `"interval"`
  (failed in (t, t_upper)).

- candset:

  column prefix for candidate set indicators, defaults to `"x"`

## Value

likelihood model object with class
`c("exp_series_md_c1_c2_c3", "series_md", "likelihood_model")`

## Details

This model satisfies the concept of a `likelihood_model` in the
`likelihood.model` package by providing the following methods:

\(1\) `loglik.exp_series_md_c1_c2_c3` (2) `score.exp_series_md_c1_c2_c3`
(3) `hess_loglik.exp_series_md_c1_c2_c3`

These are useful for doing maximum likelihood estimation, hypothesis
testing (e.g., likelihood ratio test), estimation of asymptotic sampling
distribution given data from the DGP according to the specified model,
etc.

It is designed to work well with the `likelihood_model` R package. In
particular, it is intended to be used with the `likelihood_contr_model`
object, which is a `likelihood_model` object that allows likelihood
contributions to be added for whatever data model you have in mind.

In this likelihood model, masked component data approximately satisfies
the following conditions:

C1: `Pr{K[i] in C[i]) = 1` C2:
`Pr{C[i]=c[i] | K[i]=j, T[i]=t[i]) = Pr(C[i]=c[i] | K[i]=j', T[i]=t[i])`
for any `j, j' in c[i]`. C3: masking probabilities are independent of
`theta`

As a special case, this model also includes exact component cause of
failure data where the candidate set is a singleton.

## Examples

``` r
# Create model and fit to data using generic dispatch
model <- exp_series_md_c1_c2_c3()
# solver <- fit(model)
# mle <- solver(data, par = c(1, 1, 1))
```
