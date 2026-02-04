# Random generation for a component with custom hazard/survival

Generates random samples using inverse transform sampling.

## Usage

``` r
rcomp(n, haz, surv, theta)
```

## Arguments

- n:

  number of samples to generate

- haz:

  hazard function h(t, theta)

- surv:

  survival function S(t, theta)

- theta:

  parameter vector passed to haz and surv

## Value

vector of n random samples
