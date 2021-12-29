
phi <- structure(list(
    table = list(
        series = which.min,
        parallel = which.max,
        q_out_of_m = function(q) {
            function(t) {
                #max(sort(t)[1:min(q,length(t))])
                t2 = sort(t)
                return(which.max(t[1:q]))
            }
        }),
    get_hash = function(phi) hash(phi),
    get_reverse_hash = function(h)
    {
        for (p in table)
            if (h == get_hash(p)) p
    }
))

#' system_data maps matrices of component lifetime
#' samples to data for a system with the specified
#' structure function phi.
#'
#' @param n sample size
#' @param nodes list of lifetime random variables
#'
#' @return dataframe of n observations, S x K x T, where
#'         S is system failure time, K is node failure,
#'         and T is vector of node failure times.
#' @export
#'
#' @examples
#' # generate 10 samples
#' rsystem_data(n = 10, nodes = ...)
system_data <- function(t, phi) {
    stopifnot(is.matrix(t))

    n <- nrow(t)
    m <- ncol(t)
    k <- integer(length=n)
    s <- numeric(length=n)
    for (i in 1:n)
    {
        k[i] <- phi(t[i,])
        s[i] <- t[i,k[i]]
    }

    df <- data.frame(s = s, k = k)
    df$t <- t
    attr(df,"num_nodes") <- m
    attr(df,"structure_hash") <- hash(phi)
    class(df) <- c("system_data", class(df))
    df
}


series_system_data <- function(t)
{
    system_data(t,which.min)
}

parallel_system_data <- function(t)
{
    system_data(t,which.max)
}

q_out_of_m_system_data <- function(t,q)
{
    stopifnot(ncol(t) >= q)
    stopifnot(q > 0)

    system_data(t, q_out_of_m(q,m))
}
