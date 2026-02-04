# Extract and validate model data from a masked data frame

Shared validation logic for all likelihood model methods. Checks that
the data frame is non-empty, the lifetime column exists, decodes the
candidate set matrix, and extracts the censoring indicator with
backwards compatibility.

## Usage

``` r
extract_model_data(df, lifetime, indicator, candset)
```

## Arguments

- df:

  masked data frame

- lifetime:

  column name for system lifetime

- indicator:

  column name for right-censoring indicator

- candset:

  column prefix for candidate set indicators

## Value

list with components: t (lifetimes), delta (censoring indicators), C
(candidate set matrix), m (number of components), n (number of obs)
