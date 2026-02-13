# Bernoulli candidate set model for systems with unobserved components.

Bernoulli candidate set model is a particular type of *uninformed*
model. Note that we do not generate candidate sets with this function.
See `md_cand_sampler` for that.

## Usage

``` r
md_bernoulli_cand_c1_c2_c3(
  df,
  p,
  prob = "q",
  comp = "t",
  right_censoring_indicator = "delta"
)
```

## Arguments

- df:

  masked data.

- p:

  a vector of probabilities (`p[j]` is the probability that the jth
  system will include a non-failed component in its candidate set,
  assuming the jth system is not right-censored).

- prob:

  column prefix for component probabilities, defaults to `q`, e.g.,
  `q1, q2, q3`.

- comp:

  column prefix for component lifetimes, defaults to `t`, e.g.,
  `t1, t2, t3`.

- right_censoring_indicator:

  right-censoring indicator column name. if TRUE, then the system
  lifetime is right-censored, otherwise it is observed. If NULL, then no
  right-censoring is assumed. Defaults to `delta`.

## Details

This model satisfies conditions C1, C2, and C3. The failed component
will be in the corresponding candidate set with probability 1, and the
remaining components will be in the candidate set with probability `p`
(the same probability for each component). `p` may be different for each
system, but it is assumed to be the same for each component within a
system, so `p` can be a vector such that the length of `p` is the number
of systems in the data set (with recycling if necessary).
