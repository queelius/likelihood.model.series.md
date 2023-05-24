entropy <- function(p) {
    -sum(p * log2(p))
}

hazard_general_series_helper <- function(
    hazards,
    param_counts) {

    m <- length(hazards)
    function(t, theta, log.p = FALSE) {
        i0 <- 1
        i1 <- param_counts[1]
        v <- 0
        for (j in 1:m) {
            v <- v + hazards[[j]](t,
                theta[i0:i1],
                log.p = FALSE)
            i0 <- i1 + 1
            i1 <- i1 + param_counts[j]
        }
        ifelse(log.p, log(v), v)
    }
}

survival_general_series_helper <- function(
    survivals,
    param_counts) {

    m <- length(survivals)
    function(t, theta, log.p = FALSE) {
        i0 <- 1
        i1 <- param_counts[1]
        p <- ifelse(log.p, 0, 1)
        for (j in 1:m) {
            P <- survivals[[j]](t, theta[i0:i1], log.p = log.p)
            p <- ifelse(log.p, p + P, p * P)
            i0 <- i1 + 1
            i1 <- i1 + param_counts[j]
        }
        p
    }
}
