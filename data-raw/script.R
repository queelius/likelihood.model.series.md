library(QuantPsyc)
library(MVN)

mle_exp_series <- function(md,theta0=NULL,method=mle_gradient_ascent)
{
    method(l=md_loglike_exp_series_C1_C2_C3(md),
           theta0=theta0)
}

