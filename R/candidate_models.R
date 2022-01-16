#' Decorates masked data object md with candidate sets according to
#' candidate model m0.
#'
#' @param md masked data, data frame object with column 'k' for failed component
#' and column 'w' for corresponding candidate set size.
#' @return masked data with candidate sets that model m0
#' @importFrom dplyr %>%
#' @export
md_candidate_m0 = function(md)
{
    m <- md_num_nodes(md)
    c <- tibble::as_tibble(t(apply(md,1,function(r)
    {
        x = rep(F,m)
        x[c(r["k"],sample((1:m)[-r["k"]], size=r["w"]-1, replace=F))] <- T
        x
    })))
    names(c) <- paste("c",1:m,sep=".")
    md %>% dplyr::bind_cols(c)
}

#' Decorates a masked data (input md) with candidate sets according to model m0.
#'
#' @param md masked data, a data frame object with column 'k' for failed component,
#' column 'w' for corresponding candidate set size, and column 'alpha' for
#' corresponding alpha probabilities
#' @return alpha-masked data with candidate sets that model m1
#' @importFrom dplyr %>%
#' @export
md_candidate_m1 = function(md)
{
    m <- md_num_nodes(md)
    md$test <- stats::runif(nrow(md)) < md$alpha
    c <- tibble::as_tibble(t(apply(md,1,function(r)
    {
        x = rep(F,m)
        if (r["test"])
            x[c(r["k"],sample((1:m)[-r["k"]], size=r["w"]-1, replace=F))] <- T
        else
            x[sample((1:m)[-r["k"]], size=r["w"], replace=F)] <- T
        x
    })))
    names(c) <- paste("c",1:m,sep=".")
    md %>% dplyr::bind_cols(c)
}
