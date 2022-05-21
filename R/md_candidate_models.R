#' Regular candidate model
#'
#' Decorates masked data object \code{md} with candidate sets according to
#' the regular candidate model.
#'
#' @param md masked data, data frame object with column \code{k} for failed component
#' @param m number of components in the series system
#' @return masked data
#' @importFrom dplyr %>%
#' @export
md_reg_cand = function(md,m)
{
    c <- tibble::as_tibble(t(apply(md,1,function(r)
    {
        x = rep(F,m)
        x[c(r["k"],sample((1:m)[-r["k"]], size=sample(0:(m-1),1), replace=F))] <- T
        x
    })))
    names(c) <- paste0("x",1:m)
    md %>% dplyr::bind_cols(c)
}

#' Candidate model m1
#'
#' Decorates masked data object \code{md} with candidate sets according to relaxed
#' candidate model with \code{Pr{K[i] in C[i]} = alpha}.
#'
#' @param md masked data, a data frame object with column 'k' for failed component,
#' column 'w' for corresponding candidate set size, and column 'alpha' for
#' corresponding alpha probabilities
#' @param m Integer, number of components in the series system
#' @return masked data
#' @importFrom dplyr %>%
#' @export
md_candidate_m1 = function(md,m)
{
    md$test <- stats::runif(nrow(md)) < md$alpha
    c <- tibble::as_tibble(t(apply(md,1,function(r)
    {
        x = rep(F,m)
        if (r["test"])
            x[c(r["k"],sample((1:m)[-r["k"]], size=sample(0:(m-1),1), replace=F))] <- T
        else
            x[sample((1:m)[-r["k"]], size=sample(1:m,1), replace=F)] <- T
        x
    })))
    names(c) <- paste0("x",1:m)
    md$test <- NULL
    md %>% dplyr::bind_cols(c)
}
