#' Generates a score function for an exponential series system (or
#' competing risks) with respect to parameter `theta` for masked component
#' failure (or masked competing risks) with candidate sets (risks) that
#' approximately satisfy conditions C1 and C3.
#'
#' @param md right-censored failure-time data with masked competing risks
#' @param P masking probability P{C[i] | K[i], T[i]}
#' @return score function of type `R^m -> R`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_score_exp_series_C1_C3 <- function(md, P, options, ...) {
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    defaults <- list(
        set.var = "x",
        sys.var = "t",
        delta.var = NULL)

    options <- modifyList(defaults, options)
    options <- modifyList(options, list(...))
    stopifnot(
        options$sys.var %in% colnames(md),
        is.null(options$delta.var) || options$delta.var %in% colnames(md))

    if (!is.null(options$delta.var) && options$delta.var %in% colnames(md))
        md <- md %>% filter(.data$delta == FALSE)

    t <- md[[sysvar]]
    v <- -rep(sum(t), m)

    if (!is.null(deltavar) && deltavar %in% colnames(md)) {
        # only keep the observations that were not right-censored in `md`
        md <- md %>% filter(.data$delta == FALSE)
    }

    C <- md_decode_matrix(md, setvar)
    m <- ncol(md$C)
    n <- nrow(md)

    function(theta) {
        if (length(theta) != m) stop("length(theta) != m"        
        if (any(theta <= 0)) return(NA)
        for (j in seq_len(m)) {
            for (i in seq_len(n)) {
                c <- C[i,]
                if (c[j]) {
                    denom <- 0
                    for (r in seq_len(m)) {
                        if (c[r]) {
                            denom <- denom + theta[r] * P(c, r, t[i])
                        }
                    }
                    v[j] <- v[j] + P(c, j, t[i]) / denom
                }
            }
        }
        v
    }
}






#' Generates a score function for an exponential series system (or
#' competing risks) with respect to parameter `theta` for masked component
#' failure (or masked competing risks) with candidate sets (risks) that
#' approximately satisfy conditions C1 and C3.
#'
#' @param md right-censored failure-time data with masked competing risks
#' @param P masking probability P{C[i] | K[i], T[i]}
#' @return score function of type `R^m -> R`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_score_exp_series_C1_C3 <- function(md, P, options, ...) {
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    defaults <- list(
        set.var = "x",
        sys.var = "t",
        delta.var = NULL)

    options <- modifyList(defaults, options)
    options <- modifyList(options, list(...))
    stopifnot(
        options$sys.var %in% colnames(md),
        is.null(options$delta.var) || options$delta.var %in% colnames(md))

    if (!is.null(options$delta.var) && options$delta.var %in% colnames(md))
        md <- md %>% filter(.data$delta == FALSE)

    t <- md[[sysvar]]
    v <- -rep(sum(t), m)

    if (!is.null(deltavar) && deltavar %in% colnames(md)) {
        # only keep the observations that were not right-censored in `md`
        md <- md %>% filter(.data$delta == FALSE)
    }

    C <- md_decode_matrix(md, setvar)
    m <- ncol(md$C)
    n <- nrow(md)

    function(theta) {
        if (length(theta) != m) stop("length(theta) != m"        
        if (any(theta <= 0)) return(NA)
        for (j in seq_len(m)) {
            for (i in seq_len(n)) {
                c <- C[i,]
                if (c[j]) {
                    denom <- 0
                    for (r in seq_len(m)) {
                        if (c[r]) {
                            denom <- denom + theta[r] * P(c, r, t[i])
                        }
                    }
                    v[j] <- v[j] + P(c, j, t[i]) / denom
                }
            }
        }
        v
    }
}


