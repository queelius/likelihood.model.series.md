# Quantile function for a component with custom survival function

Finds the time t such that S(t) = p using root finding. The survival
function S(t) is assumed to be monotonically decreasing from S(0) = 1 to
S(inf) = 0.

## Usage

``` r
qcomp(
  p,
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

- surv:

  survival function S(t, theta, ...)

- theta:

  parameter vector passed to surv

- t_lower:

  lower bound for search interval

- t_upper:

  upper bound for search interval (sqrt to avoid overflow)

- ...:

  additional arguments passed to surv

## Value

time t such that S(t) = p

## Examples

``` r
# Exponential survival function
surv_exp <- function(t, theta) exp(-theta * t)

# Median lifetime (50th percentile) for rate = 2
qcomp(0.5, surv = surv_exp, theta = 2.0)
#> [1] 0.3465736
```
