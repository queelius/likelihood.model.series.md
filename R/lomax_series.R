#' joint distribution of S,K,T1,...,Tm where
#'
#'     tj ~ lomax(theta_j*)
#'     s = min{T1,...,Tm}
#'     k = argmin_k {Tk : k=1,...,m}
#'
#' @param n Integer. Number of observations.
#' @param lambda Numeric vector.
#' @param kappa Numeric vector. The jth node is parameterized by theta_j := (lambda_j,kappa_j).
#' @param w Integer vector. For the ith observation, generate w_j candidates.
#' @param candidate_model the candidate model, defaults to md_candidate_m0
#' @param metadata Boolean. If TRUE writes meta-data for series system to
#'                 attributes of masked data.
#' @return masked data
#' @export
#'
#' @examples
#' md_lomax_series(n=10,lambda=c(1,2,3),kappa=c(4,5,6),w=rep(2,10))
md_lomax_series = function(n,lambda,kappa,w,candidate_model=md_candidate_m0,metadata=T)
{
    m = length(lambda)
    t <- tibble::as_tibble(matrix(
        extraDistr::rlomax(
            n*m,
            lambda=lambda,
            kappa=kappa),
        ncol=m,
        byrow=T))
    names(t) <- paste("t",1:m,sep=".")

    md <- tibble::tibble(
        s = apply(t,1,min),
        k = apply(t,1,which.min),
        w = w)
    md <- dplyr::bind_cols(md,t)
    if (!is.null(candidate_model))
        md <- candidate_model(md,m)

    if (metadata)
    {
        args <- md_func_args()
        nodes <- list()
        for (j in 1:m)
            nodes[[j]] <- list(
                "family" = "lomax",
                "index"  = c(lambda[j],kappa[j]))
        attr(md,"sim") <- list(
            "family" = toString(args[1]),
            "index" = list(lambda=lambda,kappa=kappa),
            "nodes" = nodes,
            "candidate_model" = toString(args["candidate_model"]))
        attr(md,"masked") <- c("k",paste("t",1:m,sep="."))
        attr(md,"m") <- m
        attr(md,"n") <- n
    }
    md
}
