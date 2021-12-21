library(testthat)
library(masked.data)

test_check("masked.data")


s = make_series_system(
    pdfs=c(new_univariate_pdf(function(t,x) { x*exp(-x*t) },1),
           new_univariate_pdf(function(t,x) { x*exp(-x*t) },1)),

    cdfs=c(new_univariate_cdf(function(t,x) { 1.-exp(-x*t) }),
           new_univariate_cdf(function(t,x) { 1.-exp(-x*t) })))

f1 = make.series.pdf.general(s)
f2 = make.series.pdf.general(s)
F1 = make.series.cdf.general(s)
F2 = make.series.cdf.general(s)
theta=matrix(c(1,2),nrow=2)

s2 = make.series(pdfs=c(f1,f2),
                 cdfs=c(F1,F2),
                 2)

g = make.series.pdf.general(s2)
G = make.series.cdf.general(s2)
theta2=matrix(c(1,2,3,4),nrow=2)
g(1,theta2)
G(1,theta2)

