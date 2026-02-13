# Conditional cause-of-failure probability

Returns a closure computing P(K=j \| T=t, theta) for all components j,
conditional on a specific failure time t. By Theorem 6 of the
foundational paper, this equals h_j(t; theta) / sum_l h_l(t; theta).

## Usage

``` r
conditional_cause_probability(model, ...)

# S3 method for class 'series_md'
conditional_cause_probability(model, ...)
```

## Arguments

- model:

  a likelihood model object

- ...:

  additional arguments passed to the returned closure

## Value

a function with signature `function(t, par, ...)` returning an n x m
matrix where n = length(t) and column j gives P(K=j \| T=t, theta)

## Methods (by class)

- `conditional_cause_probability(series_md)`: Default for series_md
  models via component hazard ratios (Theorem 6)
