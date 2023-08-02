


# TODOs

- test new weibull series vectorize log-like function
- use analytic score function for weibull series to MC estimate FIM (average over DGP)
- once we estimate FIM (average over DGP), we can take a very large sample
  and perform a hypothesis test to see if the MLEs are compatible
  with the asymptotic theory.


- in `generate_guo_weibull_table_2_data`, we can generate data like Table 2
  guo. It's based on that table. We estimated that, if we use
  `md_bernoulli_cand_C1_C2_C3` model with `p = 2.15`, which was the result
  of looking at the sizes of candidate sets.


- try to estimate bernoulli model in the Guo table using MLE approach. normally,
  we discard this info and don't need to estimate, but let's see what happens
  when we do. that's {C1,C2,C3}. then, try {C1,C3} for a bernoulli model where
  each probability for a candidate set can be different.