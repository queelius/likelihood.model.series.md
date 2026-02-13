# Mixture of observation schemes

Creates an observation functor that randomly selects from a set of
observation schemes for each observation. This models heterogeneous
monitoring environments where different units are observed differently.

## Usage

``` r
observe_mixture(..., weights = NULL)
```

## Arguments

- ...:

  observation functors (created by `observe_*` functions)

- weights:

  mixing probabilities (numeric vector). If NULL, uniform weights are
  used. Weights are normalized to sum to 1.

## Value

A function with signature `function(t_true)` returning a list from one
of the constituent schemes, selected randomly according to `weights`.

## Examples

``` r
obs <- observe_mixture(
  observe_right_censor(tau = 100),
  observe_left_censor(tau = 50),
  weights = c(0.7, 0.3)
)
set.seed(42)
obs(30)  # randomly selects one of the two schemes
#> $t
#> [1] 50
#> 
#> $omega
#> [1] "left"
#> 
#> $t_upper
#> [1] NA
#> 
```
