However, what if some of candidate sets are of size $w=3$ and others are, say, of size $w=2$?
    ```{r}
as_tibble(confint(md_mle_exp_series_m0(
    md_exp_series(n=n,
                  theta=theta,
                  w=c(rep(3,n/2),rep(2,n/2)),
                  candidate_model=md_candidate_m0)))) %>%
    mutate(length=.[[2]]-.[[1]])
```

We may compare this to the case when all the candidate sets are of size $w=2$ but we only have half the
sample size, $n/2$:
    ```{r}
as_tibble(confint(md_mle_exp_series_m0(
    md_exp_series(n=n/2,
                  theta=theta,
                  w=rep(2,n/2),
                  candidate_model=md_candidate_m0)))) %>%
    mutate(length=.[[2]]-.[[1]])
```

Empirically, this supports the argument that even observations that do have "candidate"
sets that include all nodes conveys information about the parameters, but if the proportion
of such complete candidate sets is too large, then the estimator is inconsistent.


#' Fisher scoring algorithm using the log-likelihood (or kernel log-like)
#' function.
#'
#' @param theta0 initial guess of theta with p components
#' @param loglike log-likelihood function of type \eqn{R^p -> R^{p \times p}}.
#' @param eps stopping condition
#' @param max_iterations
#' @parma gamma
#'
#' @return MLE estimate of theta
#' @export
md_fisher_scoring_loglike <- function(theta0,loglike,eps=1e-5,max_iterations=10000L,gamma=.1)
{
    md_fisher_scoring(
        theta0,
        function(theta) { numDeriv::hessian(function(x) { -loglike(x) },theta) },
        function(theta) { numDeriv::grad(loglike,theta) },
        eps,
        max_iterations,
        gamma)
}

#' Gradient ascent algorithm.
#'
#' @section Algorithm:
#' The algorithm is straightforward. Details here.
#'
#' @param theta0 initial guess of theta with \eqn{p} components
#' @param grad score function of type \eqn{R^p -> R^p}
#' @param eps stopping condition
#' @param max_iterations maximum number of iterations
#' @param gamma step size
#'
#' @export
md_solver_gradient_ascent <- function(theta0,grad,eps=1e-5,max_iterations=10000L,gamma=.1)
{
    n <- 1L
    repeat
    {
        theta1 <- theta0 - grad(theta0)
        if (n == max_iterations || max(abs(theta1-theta0)) < eps)
            return(list(theta.hat=theta1,iterations=n))
        theta0 <- theta1
        n <- n + 1L
    }
}







#' Kernel log-likelihood for masked data m0 for exponential series system.
#'
#' The log of the kernel of the likelihood function for masked data
#' for a series system with exponentially distributed lifetimes
#' and candidate sets that model m0.
#'
#' This is the unoptimized version, which serves as a ground-truth
#' for testing a more efficient implementation.
#'
#' @param md masked data for candidate model m0
#' @export
md_kloglike_exp_series_m0_slow <- function(md)
{
    C <- md_candidates_as_matrix(md)
    function(rate)
    {
        s <- 0.0
        for (i in 1:nrow(C))
        {
            s <- s + log(sum(rate[C[i,]]))
        }
        s - sum(md$s) * sum(rate)
    }
}








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














#Mean squared error as a function of sample size $n$:
#
#    ```{r, fig.width=6,fig.height=6,eval=F}
#Ns <- seq(1000, nrow(exp_series_data_3), 1000)
#mses.vcov <- numeric(length=length(Ns))
#mses <- numeric(length=length(Ns))
#
#for (i in 1:length(Ns))
#{
#    mses.vcov[i] <- sum(diag(vcov(md_mle_exp_series_m0(exp_series_data_3[1:Ns[i],]))))
#    mses[i] <- mean((point(md_mle_exp_series_m0(exp_series_data_3[1:Ns[i],])) - theta)^2)
#}
#
#tibble(mse=mses,mses.vcov=mses.vcov,n=Ns) %>% ggplot() +
#    geom_line(aes(x=n,y=mses,color="sample MSE")) +
#    geom_line(aes(x=n,y=mses.vcov,color="asymptotic theory")) +
#    scale_color_manual(name="Legend",values=c("sample MSE"="blue", "asymptotic theory"="red"))
#```


#However, what if some of candidate sets are of size $w=3$ and others are, say, of size $w=2$?
#    ```{r}
#as_tibble(confint(md_mle_exp_series_m0(
#    md_exp_series(n=n,
#                  theta=theta,
#                  w=c(rep(3,n/2),rep(2,n/2)),
#                  candidate_model=md_candidate_m0)))) %>%
#    mutate(length=.[[2]]-.[[1]])
#```

#We may compare this to the case when all the candidate sets are of size $w=2$ but we only have half the
#sample size, $n/2$:
#    ```{r}
#as_tibble(confint(md_mle_exp_series_m0(
#    md_exp_series(n=n/2,
#                  theta=theta,
#                  w=rep(2,n/2),
#                  candidate_model=md_candidate_m0)))) %>%
#    mutate(length=.[[2]]-.[[1]])
#```
#
#Empirically, this supports the argument that even observations that do have "candidate"
#sets that include all nodes conveys information about the parameters, but if the proportion
#of such complete candidate sets is too large, then the estimator is inconsistent.
