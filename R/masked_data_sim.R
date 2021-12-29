
#' Generate alpha-masked data from system data using probabilistic
#' model m.1. The system data must include the failed nodes.
#'
#' Model m.1 is specified in the following way. Candidate set c[j,]
#' contains the failed node and the non-failed nodes
#' are randomly selected without replacement such that
#' the cardinality of c[j,] is w[j].
#'
#' @param data a data frame with column 's' system failure time,
#'                               column 'k' failed component, and
#'                               column 't' vector of component lifetimes.
#' @param w a vector of sizes for candidate sets
#'
#' @return alpha masked data
#' @export
#'
#' @examples
#' n <- 100
#' # nodes <- list of random variables for nodes
#' # phi <- structure function, e.g., series, parallel, etc
#' data <- rsystem_data(n, nodes, phi)
#' w <- rep(2, n)
#' masked_data <- rmasked_system_data_m0(data, w)
rmasked_system_data_m0 <- function(data, w = NULL) {
    stopifnot(!is.null(data$k))
    stopifnot("num_nodes" %in% names(attributes(data)))

    k <- data$k
    m <- data$m
    n <- length(s)
    c <- matrix(nrow = n, ncol = m, F)

    if (is.null(w)) {
        w <- rep(m, n)
    }
    w <- c(w, rep(w[length(w)], n - length(w)))

    for (i in 1:n)
    {
        c[i, k[i]] <- T
        c[i, sample((1:m)[-k[i]], size = w[i] - 1, replace = F)] <- T
    }

    data$c <- c
    class(data) <- c("masked_system_data", class(data))
    attr(data,"model") <- "m0"
    data
}


#' Generate alpha-masked data from system data using probabilistic
#' model m.1. The system data must include the failed nodes.
#'
#' Model m.1 is specified in the following way. Candidate set c[j,]
#' contains failed node with probability a[j]. The non-failed nodes
#' are randomly selected without replacement. Finally,
#' sum(c[j,]) == w[j], i.e., the cardinality of the candidates are
#' pecified by vector w.
#'
#' @param data a data frame with column 's' system failure time,
#'                               column 'k' failed component, and
#'                               column 't' vector of component lifetimes.
#' @param w a vector of sizes for candidate sets
#' @param alpha a vector of alpha-probabilities
#'
#' @return alpha masked data
#' @export
#'
#' @examples
#' n <- 100
#' # nodes <- list of random variables for nodes
#' # phi <- structure function, e.g., series, parallel, etc
#' data <- rsystem_data(n, nodes, phi)
#' w <- rep(2, n)
#' alpha <- rep(.9, n)
#' masked_data <- rmasked_system_data_m1(data, alpha, w)
rmasked_system_data_m1 <- function(data, alpha = NULL, w = NULL) {
  stopifnot(!is.null(data$k))
  stopifnot("num_nodes" %in% names(attributes(data)))

  k <- data$k
  m <- data$m
  n <- length(s)
  c <- matrix(nrow = n, ncol = m, F)

  if (is.null(alpha)) {
    alpha <- rep(1, n)
  }
  alpha <- c(alpha, rep(alpha[length(alpha)], n - length(alpha)))

  if (is.null(w)) {
    w <- rep(m, n)
  }
  w <- c(w, rep(w[length(w)], n - length(w)))

  test <- stats::runif(n) < alpha
  for (i in 1:n)
  {
    if (test[i]) {
      c[i, k[i]] <- T
      c[i, sample((1:m)[-k[i]], size = w[i] - 1, replace = F)] <- T
    } else {
      c[i, sample((1:m)[-k[i]], size = w[i], replace = F)] <- T
    }
  }

  data$c <- c
  class(data) <- c("masked_system_data", class(data))
  attr(data,"model") <- "m1"

  data
}

