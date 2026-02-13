# Extract and validate model data from a masked data frame

Shared validation logic for all likelihood model methods. Checks that
the data frame is non-empty, required columns exist, decodes the
candidate set matrix, and validates observation types.

## Usage

``` r
extract_model_data(df, lifetime, omega, candset, lifetime_upper = NULL)
```

## Arguments

- df:

  masked data frame

- lifetime:

  column name for system lifetime

- omega:

  column name for observation type. Must contain character values:
  "exact", "right", "left", or "interval".

- candset:

  column prefix for candidate set indicators

- lifetime_upper:

  column name for interval upper bound (required when interval-censored
  observations are present)

## Value

list with components: t (lifetimes), omega (character vector of
observation types), C (candidate set matrix), m (number of components),
n (number of observations), t_upper (upper bounds or NULL)
