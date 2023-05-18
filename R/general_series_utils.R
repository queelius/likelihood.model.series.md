survival_general_series_helper <- function(R,nparams)
{
    m <- length(R)
    function(t,theta)
    {
        index.0 <- 1
        index.1 <- nparams[1]
        p <- 0
        for (j in 1:m)
        {
            p <- p + R[[j]](t,theta[index.0:index.1])
            index.0 <- index.1 + 1
            index.1 <- index.1 + nparams[j]
        }
        p
    }
}

hazard_general_series_helper <- function(h,nparams)
{
    m <- length(h)
    function(t,theta)
    {
        index.0 <- 1
        index.1 <- nparams[1]
        v <- 0
        for (j in 1:m)
        {
            v <- v + h[[j]](t,theta[index.0:index.1])
            index.0 <- index.1 + 1
            index.1 <- index.1 + nparams[j]
        }
        v
    }
}

generate_survival_from_hazard <- function(
    h, lowerLimit = 0, upperLimit = Inf, ...) {

    m <- length(h)
    cum_haz <- function(t) {
        hcubature(f = h, lowerLimit = lowerLimit,
            upperLimit = min(t, upperLimit), ...)$integral
    }

    function(t) exp(-cum_haz(t))
}
