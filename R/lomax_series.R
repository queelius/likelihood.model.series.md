#' rseries.lomax
#'
#' joint distribution of S,K,T1,...,Tm where
#'
#'     Tj ~ lomax(theta_j*)
#'     S = min{T1,...,Tm}
#'     K = argmin_k {Tk : k=1,...,m}
#'
#' @param n Numeric. Number of observations.
#' @param lambda Numberic vector. The j-th component has lambda=lambda_j
#' @param kappa Numberic vector. The j-th component has kappa=kappa_j
#'
#' @return matrix of n x length(lambda) component lifetimes
#' @export
#'
#' @examples
#' # generate 10 samples (10 x 3 matrix of component lifetimes)
#' t = rseries.lomax(n=10,lambda=c(1,2,3),kappa=c(4,5,6))
rseries.lomax = function(n,lambda,kappa)
{
    m = length(lambda)
    t = matrix(extraDistr::rlomax(n*m,lambda=lambda,kappa=kappa),ncol=m,byrow=T)
    K = apply(t,1,which.min)
    S = apply(t,1,min)
    data.frame(S=S,K=K,t=t)
}
