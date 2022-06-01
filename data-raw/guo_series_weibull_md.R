# The source of this data set comes from Guo, table 2, for which it is
# believed that Weibull distributed component lifetimes are a reasonable
# fit to the data and the candidate sets approximate conditions C1, C2, and C3

# The series system has m=3 components and the sample size is n=14.

guo_series_weibull_md <- read_csv("data-raw/guo_series_weibull_md.csv")
usethis::use_data(guo_series_weibull_md, overwrite=T)
