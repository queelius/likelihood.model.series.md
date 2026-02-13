# Periodic inspection observation scheme

Creates an observation functor for periodic inspections at intervals of
`delta`. Failures are bracketed between the last inspection before
failure and the first inspection after failure (interval-censored).
Systems surviving past `tau` are right-censored.

## Usage

``` r
observe_periodic(delta, tau = Inf)
```

## Arguments

- delta:

  inspection interval (positive numeric)

- tau:

  study end time (positive numeric or Inf for no right-censoring)

## Value

A function with signature `function(t_true)` returning a list with
components:

- t:

  lower bound of interval (or `tau` if right-censored)

- omega:

  "interval" or "right"

- t_upper:

  upper bound of interval (NA if right-censored)

## Examples

``` r
obs <- observe_periodic(delta = 10, tau = 100)
obs(25)   # interval: list(t = 20, omega = "interval", t_upper = 30)
#> $t
#> [1] 20
#> 
#> $omega
#> [1] "interval"
#> 
#> $t_upper
#> [1] 30
#> 
obs(150)  # right-censored: list(t = 100, omega = "right", t_upper = NA)
#> $t
#> [1] 100
#> 
#> $omega
#> [1] "right"
#> 
#> $t_upper
#> [1] NA
#> 
```
