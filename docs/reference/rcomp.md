# Random generation for a component with custom survival function

Generates random samples using inverse transform sampling.

## Usage

``` r
rcomp(n, surv, theta)
```

## Arguments

- n:

  number of samples to generate

- surv:

  survival function S(t, theta)

- theta:

  parameter vector passed to surv

## Value

vector of n random samples

## Examples

``` r
# Exponential survival function
surv_exp <- function(t, theta) exp(-theta * t)

# Generate 10 random samples with rate = 2
set.seed(123)
rcomp(10, surv = surv_exp, theta = 2.0)
#>  [1] 0.62313141 0.11893502 0.44704828 0.06220519 0.03068921 1.54440099
#>  [7] 0.31922962 0.05690974 0.29761564 0.39195764
```
