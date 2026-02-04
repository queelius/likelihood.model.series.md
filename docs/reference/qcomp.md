# Quantile function for a component with custom hazard/survival

Finds the time t such that S(t) = p using root finding. The survival
function S(t) is assumed to be monotonically decreasing from S(0) = 1 to
S(inf) = 0.

## Usage

``` r
qcomp(
  p,
  haz = NULL,
  surv,
  theta,
  t_lower = .Machine$double.eps,
  t_upper = .Machine$double.xmax^0.5,
  ...
)
```

## Arguments

- p:

  probability (quantile level), must be in (0, 1)

- haz:

  hazard function h(t, theta, ...) (currently unused but kept for API
  consistency)

- surv:

  survival function S(t, theta, ...)

- theta:

  parameter vector passed to surv

- t_lower:

  lower bound for search interval

- t_upper:

  upper bound for search interval (sqrt to avoid overflow in survival
  calculations)

- ...:

  additional arguments passed to surv

## Value

time t such that S(t) = p
