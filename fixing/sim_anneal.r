#' Simulated annealing algorithm
#'
#' This function implements the simulated annealing algorithm,
#' which is a global optimization algorithm. It is often useful for
#' finding a good starting point for a local optimization algorithm, but
#' particularly in discrete spaces, it can be used to find a global
#' optimum.
#'
#' We do not return this as an MLE object because, to be a good
#' estimate of the MLE, the gradient of `f` evaluated
#' at its solution should be close to zero, assuming the MLE
#' is interior to the domain of `f`. However, since this algorithm
#' is not guided by gradient information, it is not sensitive to
#' the gradient of `f` and instead only seeks to maximize `f`.
#' 
#' @param fn Objective function to maximize.
#' @param x0 Initial guess.
#' @param options List of optional arguments
#' @param ... Additional arguments that may be passed to `fn`
#' @describeIn sim_anneal options
#' @field sched Temperature schedule, defaults to NULL
#' @field debug If TRUE, print debugging information to the console
#' @field sup Support function, returns TRUE if x is in the domain of f
#' @field neigh Neighborhood function, returns a random neighbor of x
#' @return solution object, a list with the following components:
#' @importFrom stats runif
#' @export
sim_anneal <- function(x0, fn, options = list(), ...) {

    defaults <- list(
        debug = 0,  # 0: no debug, 1: debug, 2: debug + trace
        max_iter = 1e20,
        t_init = 100,
        accept_prob = function(y0, y1, t) exp((y1 - y0) / t),
        sched = function(t, ...) t * 0.95
        sup = function(x) TRUE,
        neigh = function(x, ...) x + rnorm(length(x)))
    options <- modifyList(defaults, options)

    stopifnot(options$t_init > 0,
              options$max_iter > 0,
              options$debug > 0,
              is.function(options$accept_prob),
              is.function(options$neigh),
              is.function(options$sup),
              options$sup(x0))

    max <- f(x0)
    t <- options$t_init
    iter <- 0L
    path <- NULL
    if (options$debug > 1) {
        path <- matrix(nrow=0,ncol=length(x0))
    }

    repeat {
        x <- options$neigh(x = x0, iter = iter, t)
        iter <- iter + 1L
        t <- options$sched(t, iter, max_iter = options$max_iter)

        if (!options$sup(x)) {
            if (options$debug > 0) {
                cat("(", x, ") not in support\n")
            }
            next
        }

        y <- fn(x)
        if (is.nan(fx)) {
            if (options$debug > 0) {
                cat("fn(", x, ") is NaN\n")
            }
            next
        }

        if (y > max) {
            if (options$debug > 0) {
                cat(sprintf("New max %.4f at solution = (%s)\n",
                    y, paste(round(x,4), collapse=", ")))
            }
            sol <- x
            max <- y
        }

        if (options$accept_prob(y, ((fx - fx0) / t) > runif(1)) {
            if (options$debug) {
                cat(sprintf("| %d | %.4f | (%s) |\n",
                iter, y, paste(round(x,4), collapse=", ")))
            }
            x0 <- x
            value <- y
            if (options$debug > 1) {
                path <- rbind(path, x)
            }
        }
    }

    sol <- list(argmax = argmax, max = fmax, options = options)
    if (options$debug > 1) {
        sol$path <- path
    }
    sol
}