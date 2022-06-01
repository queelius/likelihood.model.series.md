#' #' Construct exponential series object.
#' #'
#' #' @param rate failure rates
#' #'
#' #' @export
#' make_exp_series_dist <- function(rate)
#' {
#'     structure(list(
#'         theta=unlist(rate),
#'         num_nodes=length(rate)),
#'         class=c("exp_series_dist","series_dist","dist"))
#' }
#'
#'
#' #' Method for obtaining the variance-covariance of a \code{exp_series_dist} object.
#' #'
#' #' @param object The \code{exp_series}The object to obtain the variance of
#' #' @importFrom stats vcov
#' #' @export
#' vcov.exp_series_dist <- function(object, ...)
#' {
#'     diag(1/object$theta)^2
#' }
#'
#'
#'
#' #' Method to obtain the hazard function of
#' #' an \code{exp_series_dist} object.
#' #'
#' #' @param x The \code{exp_series_dist} object to obtain the hazard function of
#' #'
#' #' @export
#' hazard.exp_series_dist <- function(x, ...)
#' {
#'     theta <- params(x)
#'     function(t,...)
#'         ifelse(t <= 0,0,sum(theta))
#' }
#'
#'
#' #' Method to obtain the pdf of an \code{exp_series_dist} object.
#' #'
#' #' @param x The object to obtain the pdf of
#' #'
#' #' @export
#' pdf.exp_series_dist <- function(x, ...)
#' {
#'     theta <- params(x)
#'     function(t,...)
#'     {
#'         ifelse(t <= 0,0,sum(theta))
#'     }
#' }
#'
#' #' Method to sample from an \code{exp_series_dist} object.
#' #'
#' #' @param x The \code{exp_series_dist} object to sample from.
#' #' @importFrom algebraic.mle sampler
#' #' @export
#' sampler.exp_series_dist <- function(x,...)
#' {
#'     rates <- params(x)
#'     m <- length(rates)
#'     function(n=1,...)
#'         series_data(matrix(stats::rexp(m*n,rates,...),nrow=n,ncol=m))
#' }
#'
#' #' Method for obtaining the parameters of
#' #' a \code{series} distribution object.
#' #'
#' #' @param x The \code{series} object to obtain the parameters of
#' #' @importFrom algebraic.mle params
#' #' @export
#' params.series <- function(x, ...)
#' {
#'     x$theta
#' }




#' #' @param md masked data
#' #' @param h list of hazard functions
#' #' @param measure a measure function
#' #' @importFrom md.tools md_decode_matrix
#' #' @importFrom dplyr %>%
#' #' @importFrom dplyr bind_cols
#' #' @importFrom stats runif
#' #' @export
#' md_filter_measure <- function(md,h,measure)
#' {
#'     md_decode_matrix(md, "pk")
#'     h.series <- hazard_general_series(h)
#'     t <- md$t
#'     n <- nrow(md)
#'
#'     pk <- matrix(rep(NA),m*n,ncol=m)
#'     right_censoring <- "delta" %in% colnames(md)
#'     for (i in 1:n)
#'     {
#'         if (!right_censoring || (right_censoring && !md$delta[i]))
#'         {
#'             for (j in 1:m)
#'                 pk[i,j] <- h[[j]](t[i]) / h.C(t[i])
#'         }
#'     }
#'     md %>% add_column(
#'         }
#' #' Construct exponential series object.
#' #'
#' #' @param rate failure rates
#' make_exp_series_dist <- function(rates)
#' {
#'     structure(list(
#'         theta=unlist(rates),
#'         num_comp=length(rates)),
#'         class=c("exp_series_dist","series_dist","exp_dist","dist"))
#' }
#'
#' #' Method for obtaining the variance-covariance of a \code{exp_series_dist} object.
#' #'
#' #' @param object The \code{exp_series}The object to obtain the variance of
#' #' @importFrom stats vcov
#' vcov.exp_series_dist <- function(object, ...)
#' {
#'     1/sum(params(x))^2
#' }
#'
#' #' Method to obtain the hazard function of
#' #' an \code{exp_series_dist} object.
#' #'
#' #' @param x The \code{exp_series_dist} object to obtain the hazard function of
#' hazard.exp_series_dist <- function(x, ...)
#' {
#'     rate <- sum(params(x))
#'     function(t,...)
#'         ifelse(t < 0,0,rate)
#' }
#'
#'
#' #' Method to obtain the pdf of an \code{exp_series_dist} object.
#' #'
#' #' @param x The object to obtain the pdf of
#' pdf.exp_series_dist <- function(x, ...)
#' {
#'     rate <- sum(params(x))
#'     function(t,...)
#'     {
#'         ?dexp()
#'     }
#' }
#'
#' #' Method to sample from an \code{exp_series_dist} object.
#' #'
#' #' @param x The \code{exp_series_dist} object to sample from.
#' #' @importFrom algebraic.mle sampler
#' sampler.exp_series_dist <- function(x,...)
#' {
#'     rates <- params(x)
#'     m <- length(rates)
#'     function(n=1,...)
#'         stats::rexp(n,sum(rates),...)
#' }
#'
#' #' Method for obtaining the parameters of
#' #' a \code{series} distribution object.
#' #'
#' #' @param x The \code{series} object to obtain the parameters of
#' #' @importFrom algebraic.mle params
#' params.series <- function(x, ...)
#' {
#'     x$theta
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' md_loglike_weibull_series2 <- function(md,cand="x")
#' {
#'     n <- nrow(md)
#'     ud <- NULL
#'     if (delta %in% colnames(md))
#'         ud <- filter(md,delta==F)
#'
#'     un <- length(ud)
#'     C <- md_decode_matrix(ud,cand)
#'     m <- ncol(C)
#'     stopifnot(m > 0)
#'
#'     function(theta)
#'     {
#'         # theta should be a parameter vector of length 2*m
#'         # todo: make first m values of `theta` be for scale parameters
#'         #            second m values (m+1,...,2m) be for shape parameters
#'         scales <- theta[(0:(m-1)*2)+1]
#'         shapes <- theta[(1:m)*2]
#'         s <- 0
#'
#'         for (i in 1:n)
#'         {
#'             acc <- 0
#'             for (j in 1:m)
#'                 acc <- acc + (md$s[i]/scales[j])^shapes[j]
#'             s <- s + log(acc)
#'         }
#'
#'         #function(r)
#'         #{
#'         #    sum(log((r["s"]/scales)^shapes))
#'         #}
#'
#'
#'         for (i in 1:un)
#'         {
#'             acc <- 0
#'             c <- (1:m)[C[i,]]
#'             for (j in c)
#'                 acc <- acc + shapes[j]/scales[j]*(ud$s[i]/scales[j])^(shapes[j]-1)
#'             s <- s + log(acc)
#'         }
#'
#'         s
#'     }
#' }
#'
#' #' Constructs a pdf (pmf) object for the conditional probability of component
#' #' failure \code{Pr{K[i]=j|C[i]=c[i],T[i]=t[i]) = h_j(t[i])/h(t[i]) I(j in c[i])}.
#' #' in an exponential series system for the regular candidate model.
#' #'
#' #' @param theta parameter value of \code{exp_series_dist}
#' #' @export
#' md_exp_series_component_failure_decorator_C1 <- function(md,theta)
#' {
#'     if (is.matrix(theta))
#'         theta <- as.vector(theta)
#'     theta <- unlist(theta)
#'     n <- nrow(md)
#'     m <- length(theta)
#'     md$C <- md.tools::md_decode_matrix(md,"x")
#'
#'     K <- matrix(rep(NA,m*n),ncol=m)
#'     for (i in 1:n)
#'     {
#'         for (k in 1:m)
#'         {
#'             if (md$s[i] <= 0 || !(k %in% (1:m)[md$C[i,]]))
#'                 K[i,k] <- 0
#'             else
#'                 K[i,k] <- theta[k]/sum(theta[md$C[i,]])
#'         }
#'     }
#'
#'     K <- tibble::as_tibble(K)
#'     colnames(K) <- paste0("k.",1:m)
#'     md <- dplyr::bind_cols(md,K)
#'     md$C <- NULL
#'     md
#' }




#'
#'
#'
#'
#'
#'
#' #' Compute the covariance matrix from the given masked data estimate.
#' #'
#' #' Sampling distribution of the MLE is a multivariate normal with mean
#' #' given by the true parameter value and, asymptotically, a covariance
#' #' given by the inverse of the Fisher information matrix.
#' #'
#' #' @param object The variance-covariance matrix of the estimator to obtain
#' #' @param ... Additional arguments to pass.
#' #' @importFrom stats vcov
#' #' @export
#' vcov.md_mle <- function(object, ...)
#' {
#'     object$sigma
#' }
#'
#' #' Method to obtain the confidence intervals of the parameter values of a
#' #' masked data estimator, \code{md_estimate}.
#' #'
#' #' @param object The \code{md_estimate} object to compute the confidence intervals for
#' #' @param parm Unused
#' #' @param level Confidence level, defaults to 0.95 (alpha=.05)
#' #' @param ... Additional arguments to pass.
#' #' @importFrom stats confint
#' #' @export
#' confint.md_mle <- function(object, parm=NULL, level=0.95, ...)
#' {
#'     V <- vcov(object)
#'     theta.hat <- point(object)
#'     p <- length(theta.hat)
#'     q <- stats::qnorm(level)
#'
#'     ci <- matrix(rep(NA,p*2),p,2)
#'     colnames(ci) <- c(paste((1-level)/2*100,"%",sep=""),
#'                       paste((1-(1-level)/2)*100,"%",sep=""))
#'     for (j in 1:p)
#'     {
#'         ci[j,1] <- theta.hat[j] - q * sqrt(V[1,1])
#'         ci[j,2] <- theta.hat[j] + q * sqrt(V[1,1])
#'     }
#'     ci
#' }
#'
#' #' Method to obtain the point estimate of
#' #' a masked data estimator, \code{md_mle}.
#' #'
#' #' @param x The \code{md_mle} object to obtain the MLE point estimate of
#' #' @param ... Additional arguments to pass.
#' #' @importFrom algebraic.mle point
#' #' @export
#' point.md_mle <- function(x, ...)
#' {
#'     x$theta.hat
#' }
#'
#' #' Method to obtain the fisher information matrix of an \code{md_mle}.
#' #'
#' #' @param x The \code{md_mle} object to obtain the fisher information of
#' #' @param ... Additional arguments to pass.
#' #'
#' #' @importFrom algebraic.mle fisher_info
#' #' @export
#' fisher_info.md_mle <- function(x, ...)
#' {
#'     x$info
#' }
#'
#' #' Method to obtain the sampler for an \code{md_mle} object.
#' #'
#' #' @param x The \code{md_mle} object to create a sampling procedure from
#' #' @param ... Additional arguments to pass.
#' #' @importFrom algebraic.mle sampler
#' #' @importFrom mvtnorm rmvnorm
#' #' @export
#' sampler.md_mle <- function(x, ...)
#' {
#'     function(n) rmvnorm(n,point(x),vcov(x))
#' }
