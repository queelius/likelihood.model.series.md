#' Bernoulli candidate model that is a function of component
#' cause of failure.
#'
#' Bernoulli candidate model that is a function of component
#' cause of failure \code{k} and parameters \code{m} and \code{p},
#' where if component cause of failure is indexed by \code{j}, then \code{j} is
#' in the candidate set and otherwise the component indexed by \code{j} is in
#' the candidate set with probability specified by \code{p}, which may either
#' by a procedure for sampling from a distribution with a codomain \code{[0,1]}
#' or a constant function.
#'
#' @param md masked data, has a column \code{k} for the failed component index.
#' @param m integer, number of components in the candidate set.
#' @param p a function or procedure \code{p : integer -> [0,1]}, defaults to
#'          \code{function(n) runif(n)}.
#' @importFrom md.tools md_decode_matrix
#' @importFrom dplyr %>%
#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#' @importFrom stats runif
#' @export
md_bernoulli_candidate_C1_C2_C3 <- function(md,m,p=function(n) runif(n))
{
    stopifnot(!is.null(md$k))
    n <- nrow(md)
    stopifnot(n > 0)

    x <- matrix(NA,nrow=n,ncol=m)
    u <- matrix(runif(m*n),nrow=n)
    gam <- p(n)

    for (i in 1:n)
    {
        for (j in 1:m)
        {
            x[i,j] <- ifelse(md$k[i]==j,
                             T,
                             u[i,j] < gam[i])
        }
    }

    x <- tibble::as_tibble(x)
    colnames(x) <- paste0("x",1:m)
    md %>% dplyr::bind_cols(x)
}

md_block_candidate_m3 <- function(md)
{
    block <- function(k)
    {
        if (k == 1)
            return(c(T,T,F))
        if (k == 2)
            return(c(T,T,F))
        if (k == 3)
        {
            if (runif(1) < 0.1)
                return(c(T,T,T))
            else
                return(c(F,F,T))
        }
    }

    for (i in 1:n)
        x[i,] <- block(md$k[i])

    x <- tibble::as_tibble(x)
    colnames(x) <- paste0("x",1:3)
    md %>% dplyr::bind_cols(x)
}
