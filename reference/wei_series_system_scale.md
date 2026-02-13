# System scale parameter for homogeneous Weibull series

For a series system with Weibull components sharing shape k but with
individual scales, the system lifetime is itself Weibull with shape k
and this computed scale.

## Usage

``` r
wei_series_system_scale(k, scales)
```

## Arguments

- k:

  common shape parameter

- scales:

  vector of component scale parameters

## Value

system scale parameter

## Examples

``` r
# 3-component system with common shape 1.2
wei_series_system_scale(k = 1.2, scales = c(1000, 900, 850))
#> [1] 365.1326
```
