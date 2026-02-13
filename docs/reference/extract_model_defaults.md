# Extract model column name defaults

Helper function to extract default column names from a likelihood model
object. Used by all model methods to avoid repeating the same pattern.

## Usage

``` r
extract_model_defaults(model)
```

## Arguments

- model:

  likelihood model object with lifetime, lifetime_upper, omega, candset
  fields

## Value

list with lifetime, lifetime_upper, omega, candset defaults
