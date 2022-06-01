# Retrieve the function arguments.
md_func_args <- function(...)
{
    call <- evalq(match.call(expand.dots = F), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))

    for(i in setdiff(names(formals), names(call)))
        call[i] <- list( formals[[i]] )

    match.call(sys.function(sys.parent()), call)
}

# mean absolute deviation
mad <- function(x)
{
    u <- mean(x)
    n <- length(x)
    sum(abs(x-u)) / n
}

# linearly rescale numbers in x each in the range [x_low,x_high] to an output
# in the range [out_low, out_high].
rescale <- function(x,x_low=0,x_high=1,y_low=0,y_high=1)
{
    y_low + ((y_high - y_low) / (x_high - x_low)) * (x - x_low)
}

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

