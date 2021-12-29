#' rseries_system_data generates the joint distribution
#' of a series system with the indicated nodes.
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
#' rseries_system_data(n = 10, nodes = ...)
rseries_system_data <- function(n, nodes) {
    df <- rsystem_data(n,nodes,which.min)
    class(df) = c("series_system_data", class(df))
    df
}

#' rpar_system_data generates data for the joint distribution
#' of a parallel system with the indicated nodes.
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
#' rseries_system_data_from_nodes(n = 10, nodes = ...)
rseries_system_data_from_nodes <- function(n, nodes) {
    df <- rsystem_data_from_nodes(n,nodes,which.max)
    class(df) = c("parallel_system_data", class(df))
    df
}

#' rsystem_data_from_nodes generates/simulates
#' data from a system with the indicated nodes
#' and structure function phi.
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
#' rsystem_data_from_nodes(n = 10, nodes = ...)
rsystem_data_from_nodes <- function(n, nodes, phi) {
    t <- NULL
    for (node in nodes) {
        t <- cbind(t, node(n))
    }
    system_data(t,phi)
}


