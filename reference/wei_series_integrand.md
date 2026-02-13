# Integrand for numerical integration of Weibull series likelihood

Computes h_c(t) \* S(t) where h_c = sum_j in c h_j(t) and S(t) =
exp(-sum_l H_l(t)). Used for left-censored and interval-censored
observations in the heterogeneous Weibull model.

## Usage

``` r
wei_series_integrand(t, shapes, scales, cind)
```

## Arguments

- t:

  time values (vector, for use with stats::integrate)

- shapes:

  shape parameters for all components

- scales:

  scale parameters for all components

- cind:

  logical vector indicating which components are in candidate set

## Value

vector of integrand values
