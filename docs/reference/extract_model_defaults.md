# Extract model column name defaults

Helper function to extract default column names from a likelihood model
object. Used by all model methods to avoid repeating the same pattern.

## Usage

``` r
extract_model_defaults(model)
```

## Arguments

- model:

  likelihood model object with lifetime, indicator, candset fields

## Value

list with lifetime, indicator, candset defaults
