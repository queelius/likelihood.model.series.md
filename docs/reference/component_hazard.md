# Component hazard function

Returns a closure computing the hazard function h_j(t; theta) for the
j-th component. The returned function takes the full parameter vector
`par` and extracts the relevant component parameters internally.

## Usage

``` r
component_hazard(model, j, ...)
```

## Arguments

- model:

  a likelihood model object

- j:

  component index (integer from 1 to m)

- ...:

  additional arguments passed to the returned closure (e.g., covariates
  for proportional hazards extensions)

## Value

a function with signature `function(t, par, ...)` computing h_j(t)
