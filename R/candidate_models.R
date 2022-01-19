#' Candidate model m0
#'
#' Decorates masked data object \code{md} with candidate sets according to candidate
#' model m0.
#'
#' Specifically, the candidate sets are generated according to the alpha-masked
#' model, where C_i contains k_i and w_i-1 nodes randomly selected without
#' replacement from {1,...,m} - {k_i}.
#'
#' @param md masked data, data frame object with column 'k' for failed component
#' and column 'w' for corresponding candidate set size.
#' @param m Integer, number of nodes in the series system
#' @return masked data with candidate sets that model m0
#' @importFrom dplyr %>%
#' @export
md_candidate_m0 = function(md,m)
{
    c <- tibble::as_tibble(t(apply(md,1,function(r)
    {
        x = rep(F,m)
        x[c(r["k"],sample((1:m)[-r["k"]], size=r["w"]-1, replace=F))] <- T
        x
    })))
    names(c) <- paste("c",1:m,sep=".")
    md %>% dplyr::bind_cols(c)
}

#' Candidate model m1
#'
#' Decorates masked data object \code{md} with candidate sets according to candidate
#' model m1.
#'
#' Specifically, the candidate sets are generated according to the alpha-masked
#' model, where with probability alpha_i, C_i contains k_i and w_i-1 nodes
#' randomly selected without replacement from {1,...,m} - {k_i} and with
#' probability 1-alpha_i, C_i contains w_i nodes randomly selected without
#' replacement from {1,...,m} \ { k_i }.
#'
#' @param md masked data, a data frame object with column 'k' for failed component,
#' column 'w' for corresponding candidate set size, and column 'alpha' for
#' corresponding alpha probabilities
#' @param m Integer, number of nodes in the series system
#' @param with_seed Integer, seeds PRNG to facilitate reproducible data/research
#' @return alpha-masked data with candidate sets that model m1
#' @importFrom dplyr %>%
#' @export
md_candidate_m1 = function(md,m)
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
    md$test <- NULL
    md %>% dplyr::bind_cols(c)
}
