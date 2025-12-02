# Sample candidate sets for systems with unobserved components.

Candidate set generator. Requires columns for component probabilities
e.g., `q1,...,qm` where `qj` is the probability that the jth component
will be in the corresponding candidate set generated for that
observation in the `md` table.

## Usage

``` r
md_cand_sampler(df, prob = "q", candset = "x")
```

## Arguments

- df:

  (masked) data frame

- prob:

  column prefix for component probabilities, defaults to `q`, e.g.,
  `q1, q2, q3`.

- candset:

  column prefix for candidate sets (as Boolean matrix), defaults to `x`,
  e.g., `x1, x2, x3`.
