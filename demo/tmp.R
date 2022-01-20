#' Construct a new system
#'
#' Constructs a new lifetime system with nodes that
#' have the specified set of log-densities.
#'
#' @param log_pdfs parametric log-densities for components,
#'             make sure that each log-density accepts
#'             the same number parameter vector.
#'
#' @return Object representing the system without any
#'         specified structure.
#' @export
new_system = function(log_pdfs)
{
    #stopifnot(is_pdf(pdfs))
    #stopifnot(is_cdf(cdfs))

    pdfs[[1]]$param.dim

    structure(list(pdfs=pdfs,cdfs=cdfs),
              param.dim=c(length(pdfs),p),
              class=c("series_system","random_variable"))
}


















#' bootstrap of the mle's sampling distribution
#'
#' @param masked.data sample of masked data
#' @param r replicates
#'
#' @return mle and bootstrap estimate of its covariance matrix
#' @export
# rseries.exp.mle.cov.bootstrap = function(masked.data,r)
# {
#     n = nrow(masked.data)
#     rate.mle = series.exp.mle(masked.data)
#     rate.bs = NULL
#     for (i in 1:r)
#     {
#         d = rmasked.data(n,rate.mle)
#         rate.bs = rbind(rate.bs,series.exp.mle(d))
#     }
#     list(mle=rate.mle,Sigma=cov(rate.bs))
# }
#
#
#
#
#
#
# mse = function(theta,n) { length(theta)*sum(theta)/n }
#
#
# info <- function(l,n)
# {
#     l1 = l[1]
#     l2 = l[2]
#     l3 = l[3]
#
#     A = matrix(c(1/(l1+l2) + 1/(l1+l3), 1/(l1+l2), 1/(l1+l3),
#                  1/(l1+l2), 1/(l1+l2) + 1/(l2+l3), 1/(l2+l3),
#                  1/(l1+l3), 1/(l2+l3), 1/(l1+l3) + 1/(l2+l3)),byrow=T,nrow=3)
#     n/(2.0*(l1+l2+l3))*A
# }
#
#
# #delta=theta.mle - theta.star
#
# mse.mle = sum(diag(cov.mle))
# mse.a = mse(theta.mle,n)
# mse.star = mse(theta.star,n)
#
# #### sufficient statistics
#



#' Construct a new series system. We do not make
#' any assumptions about the parametric family
#' of the nodes in series.
#'
#' @description
#' Constructs a new series system from the specified
#' list of univariate distribution functions, the
#' pdfs and cdfs. The distribution functions can be
#' a mix of any parametric families.
#'
#' @param pdfs Parametric pdf functions for components.
#'             each should accept a failure time and a
#'             parameter matrix theta, which denotes
#'                 fk(t|theta).
#'
#' @param cdfs Parametric cdf functions for components.
#'             each should accept a failure time and a
#'             parameter matrix theta, which denotes
#'                 Fk(t|theta).
#'
#' @return Object representing the series system.
#' @export
new_series_system = function(pdfs,cdfs)
{
    #stopifnot(is_pdf(pdfs))
    #stopifnot(is_cdf(cdfs))

    pdfs[[1]]$param.dim

    structure(list(pdfs=pdfs,cdfs=cdfs),
              param.dim=c(length(pdfs),p),
              class=c("series_system","random_variable"))
}

new_series_system_pdf <- function(series) {
    # stopifnot(is_series_system(series))

    m <- attributes(series)$param.dim[1]
    p <- attributes(series)$param.dim[2]

    structure(
        function(t, theta) {
            # stopifnot(is.matrix(theta))
            # stopifnot(nrow(theta) == m)
            # stopifnot(ncol(theta) == p)
            theta <- as.matrix(theta, nrow = m, ncol = p)
            d <- 0.
            for (j in 1:m)
            {
                dj <- series$pdfs[[j]](t, theta[j, ])
                for (k in (1:m)[-j]) {
                    dj <- dj * (1 - series$cdfs[[k]](t, theta[k, ]))
                }
                d <- d + dj
            }
            d
        },
        param.dim = c(m, p),
        class = c("series_system_pdf", "pdf", "function")
    )
}

new_series_system_cdf = function(series)
{
    # stopifnot(is.series.system(series))

    m = attributes(series)$param.dim[1]
    p = attributes(series)$param.dim[2]

    structure(
        function(t,theta)
        {
            #stopifnot(is.matrix(theta))
            #stopifnot(nrow(theta) == m)
            #stopifnot(ncol(theta) == p)
            theta=as.matrix(theta,nrow=m,ncol=p)
            d = 1.
            for (j in 1:m)
                d = d * (1-series$cdfs[[j]](t,theta[j,]))
            return(1-d)
        },
        param.dim=c(m,p),
        class=c("series_system_cdf","cdf","function")
    )
}






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

















#' Compute the covariance matrix from a sample
#' of MLEs.
#'
#' The asymptotic sampling distribution of the MLE
#' is a multivariate normal with mean
#' given by the true parameter value
#' and a covariance given by the inverse of its
#' Fisher information matrix.
#'
#' However, if we are not at the asymptotic
#' limit, we may estimate the sampling distribution
#' if we have ...
#'
#' @param mles vector of MLEs
#' @export
md_mle_covariance <- function(mles)
{
    NA
}






lam <- c(3,4,5)
h <- list(
    function(t) { lam[1] },
    function(t) { lam[2] },
    function(t) { lam[3] })

fk <- md_series_node_failure_m0(h)




est <- md_mle_exp_series_m0(exp_series_data_1)
mle <- point(est)
md_loglik <- md_kloglike_lomax_series_m0_ref(lomax_series_data_1)
md_l <- function(x,y) { md_loglik(c(x,y,mle[3])) }

# define the drawing grid
# create all possible combinations of x and y
inputs <- dplyr::tibble(x=seq(mle[1]-1,mle[1]+1,0.05),y=seq(mle[2]-1,mle[2]+1,0.05)) %>% purrr::cross_df()

n <- nrow(inputs)
surf <- inputs
surf$z <- numeric(n)
for (i in 1:n)
{
    surf$z[i] <- md_l(surf$x[i],surf$y[i])
}

surf %>% ggplot(aes(x=x,y=y,z=z)) + geom_contour()

