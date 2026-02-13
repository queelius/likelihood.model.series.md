# Left-censoring observation scheme (single inspection)

Creates an observation functor for a single-inspection design. If the
system has already failed by inspection time `tau`, we know it failed
before `tau` but not exactly when (left-censored). If it is still
running, we know it survived past `tau` (right-censored).

## Usage

``` r
observe_left_censor(tau)
```

## Arguments

- tau:

  inspection time (positive numeric)

## Value

A function with signature `function(t_true)` returning a list with
components:

- t:

  inspection time `tau`

- omega:

  "left" if failed before `tau`, "right" otherwise

- t_upper:

  NA (not used for this scheme)

## Examples

``` r
obs <- observe_left_censor(tau = 100)
obs(50)   # left-censored: list(t = 100, omega = "left", t_upper = NA)
#> $t
#> [1] 100
#> 
#> $omega
#> [1] "left"
#> 
#> $t_upper
#> [1] NA
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
