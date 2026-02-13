# Marginal cause-of-failure probability

Returns a closure computing P(K=j \| theta) for all components j,
marginalized over the system failure time T. By Theorem 5 of the
foundational paper, this equals E_T\[P(K=j \| T, theta)\].

## Usage

``` r
cause_probability(model, ...)

# S3 method for class 'series_md'
cause_probability(model, ...)
```

## Arguments

- model:

  a likelihood model object

- ...:

  additional arguments passed to the returned closure

## Value

a function with signature `function(par, ...)` returning an m-vector
where element j gives P(K=j \| theta)

## Details

The default method uses Monte Carlo integration via
[`rdata()`](https://queelius.github.io/likelihood.model/reference/rdata.html).

## Methods (by class)

- `cause_probability(series_md)`: Default for series_md models via Monte
  Carlo integration (Theorem 5)
