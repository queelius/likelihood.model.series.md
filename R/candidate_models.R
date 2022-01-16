#' rcandidates.m0 : m0 model (1-masked data)
#'
#' @param md a data frame with column 'k' for failed component and
#'           column 'w' for corresponding candidate set size.
#' @return masked data with candidate sets that model m0
#' @export
#'
#' @examples
#' data = rcandidates.m0(md)
md_candidates_m0 = function(md)
{
    c <- as_tibble(t(apply(md,1,function(r)
    {
        x = rep(F,m)
        x[c(r["k"],sample((1:m)[-r["k"]], size=r["w"]-1, replace=F))] <- T
        x
    })))
    names(c) <- paste("c",1:m,sep=".")
    md %>% bind_cols(c)
}

#' Generate masked data from series system data
#'
#' @param md a data frame with column 'k' for failed component,
#'           column 'w' for corresponding candidate set size,
#'           and column 'alpha' for corresponding alpha probabilities
#'
#' @return alpha-masked data with candidate sets that model m1
#' @export
#'
#' @examples
#' md <- rcandidates.m1(md)
md_candidates_m1 = function(md)
{
    md$test <- stats::runif(n) < md$alpha
    c <- as_tibble(t(apply(md,1,function(r)
    {
        x = rep(F,m)
        if (r["test"])
            x[c(r["k"],sample((1:m)[-r["k"]], size=r["w"]-1, replace=F))] <- T
        else
            x[sample((1:m)[-r["k"]], size=r["w"], replace=F)] <- T
        x
    })))
    names(c) <- paste("c",1:m,sep=".")
    md %>% bind_cols(c)
}
