# Mark columns as latent in a masked data frame

Sets the `"latent"` attribute on a data frame, recording which columns
represent unobserved (latent) variables.

## Usage

``` r
md_mark_latent(md, vars)
```

## Arguments

- md:

  data frame to modify

- vars:

  character vector of column names to mark as latent

## Value

the data frame with updated latent attribute
