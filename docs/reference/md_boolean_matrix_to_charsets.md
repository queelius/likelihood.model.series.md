# Convert Boolean candidate set columns to character set notation

Replaces Boolean matrix columns (e.g., `x1, x2, x3`) with a single
character column showing set notation like `{1, 3}`.

## Usage

``` r
md_boolean_matrix_to_charsets(df, setvar = "x", cname = NULL, drop_set = FALSE)
```

## Arguments

- df:

  data frame containing Boolean matrix columns

- setvar:

  column prefix for the Boolean matrix (default `"x"`)

- cname:

  name for the new character column (default: same as `setvar`)

- drop_set:

  if TRUE, remove the original Boolean columns (default FALSE)

## Value

data frame with character set column added
