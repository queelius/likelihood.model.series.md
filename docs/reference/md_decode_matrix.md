# Decode a matrix from prefixed columns in a data frame

Extracts columns matching the pattern `var1, var2, ...` or
`var.1, var.2, ...` from a data frame and returns them as a numeric
matrix, ordered by index.

## Usage

``` r
md_decode_matrix(df, var)
```

## Arguments

- df:

  data frame containing the matrix columns

- var:

  character prefix for the column names

## Value

a matrix, or NULL if no matching columns are found
