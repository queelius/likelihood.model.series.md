# Right-censoring observation scheme

Creates an observation functor that applies right-censoring at time
`tau`. Systems that fail before `tau` are observed exactly; systems
surviving past `tau` are right-censored.

## Usage

``` r
observe_right_censor(tau)
```

## Arguments

- tau:

  censoring time (positive numeric)

## Value

A function with signature `function(t_true)` returning a list with
components:

- t:

  observed time

- omega:

  "exact" or "right"

- t_upper:

  NA (not used for this scheme)

## Examples

``` r
obs <- observe_right_censor(tau = 100)
obs(50)   # exact: list(t = 50, omega = "exact", t_upper = NA)
#> $t
#> [1] 50
#> 
#> $omega
#> [1] "exact"
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
