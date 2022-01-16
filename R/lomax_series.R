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
md_lomax_series_m0 = function(n,lambda,kappa,w)
{
    m = length(lambda)
    t <- as_tibble(matrix(extraDistr::rlomax(n*m,lambda=lambda,kappa=kappa), ncol=m, byrow=T))
    names(t) <- paste("t",1:m,sep=".")

    md <- tibble(s = apply(t,1,min),
                 k = apply(t,1,which.min),
                 w = w)
    md <- md_candidates_m0(md)
    md <- bind_cols(md,t)

    # set up attribute properties
    nodes <- list()
    for (j in 1:m)
        nodes[[j]] <- list("family" = "lomax",
                           "index"  = c(lambda[j],kappa[j]))

    attr(md,"sim") <- list("theta" = list(lambda=lambda,kappa=kapp),
                           "nodes" = nodes,
                           "model" = "m0")
    attr(md,"masked") <- c("k",paste("t",1:m,sep="."))

    md
}
