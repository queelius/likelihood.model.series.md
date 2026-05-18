#' Exponential series systems (legacy convenience API)
#'
#' These functions provide a base-R-style distribution API
#' (`d`/`p`/`q`/`r`/hazard/surv/mean) for exponential series systems and
#' are preserved for backwards compatibility. New code should construct
#' a [dist.structure::exp_series()] object, which provides the same dist
#' interface plus broader interop (importance measures, structural
#' importance, system signature, dual, compositional operations).
#'
#' A series system of m independent exponential components with rates
#' \eqn{\lambda_1, \ldots, \lambda_m} has an exponential system lifetime
#' with rate \eqn{\lambda = \sum \lambda_j}.
#'
#' @name Exponential-series-system
#' @keywords internal
NULL

#' Quantile function for exponential series.
#'
#' Thin wrapper around [stats::qexp()] with `rate = sum(rates)`. Equivalent
#' to `algebraic.dist::inv_cdf(dist.structure::exp_series(rates))(p)`.
#'
#' @param p Vector of quantiles.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param lower.tail Logical, if TRUE (default), probabilities are P(X <= x),
#'                   otherwise, P(X > x).
#' @param log.p Logical, if TRUE, vector of probabilities `p` are given as `log(p)`.
#' @return Quantiles corresponding to the given probabilities `p`.
#' @note Preserved for backwards compatibility. New code: use
#'   [dist.structure::exp_series()] then [algebraic.dist::inv_cdf()].
#' @importFrom stats qexp
#' @export
#' @examples
#' rates <- c(0.5, 0.3, 0.2)
#' qexp_series(0.5, rates)
#' qexp_series(c(0.25, 0.5, 0.75), rates)
qexp_series <- function(p, rates, lower.tail = TRUE, log.p = FALSE) {
    # Series of m independent exponentials is Exp(sum(rates)); the
    # closed-form quantile is qexp at the aggregate rate. The
    # dist.structure::exp_series(rates) aggregate rate is the canonical
    # source.
    qexp(p, dist.structure::exp_series(rates)$total_rate, lower.tail, log.p)
}

#' Random number generation for exponential series.
#'
#' Generates random variates from an exponential series distribution.
#' The default path (`keep_latent = FALSE`) is equivalent to
#' `algebraic.dist::sampler(dist.structure::exp_series(rates))(n)`. The
#' `keep_latent = TRUE` path additionally returns the latent component
#' lifetimes (sampled independently from each component).
#'
#' @param n Integer, number of observations to generate.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param keep_latent Logical; if TRUE, returns a matrix with system lifetimes
#'                    in the first column and individual component lifetimes in
#'                    subsequent columns. If FALSE (default), returns only
#'                    system lifetimes.
#' @return If `keep_latent = FALSE`, a vector of random variates from the
#'         exponential series distribution. If `keep_latent = TRUE`, a matrix
#'         with system lifetime in the first column and component lifetimes
#'         in columns 2 through m+1.
#' @note Preserved for backwards compatibility. New code: use
#'   [dist.structure::exp_series()] with [algebraic.dist::sampler()].
#' @importFrom stats rexp
#' @export
#' @examples
#' set.seed(123)
#' rexp_series(5, rates = c(0.5, 0.3, 0.2))
#' rexp_series(3, rates = c(0.5, 0.3, 0.2), keep_latent = TRUE)
rexp_series <- function(n, rates, keep_latent = FALSE) {
    if (!keep_latent) {
        return(rexp(n, dist.structure::exp_series(rates)$total_rate))
    }
    # Latent component lifetimes path: sample each component separately.
    m <- length(rates)
    data <- matrix(NA, nrow = n, ncol = m)
    for (j in seq_len(m)) {
        data[, j] <- rexp(n, rates[j])
    }
    cbind(apply(data, 1, min), data)
}

#' Density function for exponential series.
#'
#' Thin wrapper. Equivalent to
#' `algebraic.dist::density(dist.structure::exp_series(rates))(t, log = log)`.
#'
#' @param t series system lifetime
#' @param rates rate parameters for exponential component lifetimes
#' @param log return the log of the pdf
#' @return Density values for the exponential series distribution.
#' @note Preserved for backwards compatibility. New code: use
#'   [dist.structure::exp_series()] then `algebraic.dist::density()`.
#' @importFrom stats dexp
#' @export
#' @examples
#' rates <- c(0.5, 0.3, 0.2)
#' dexp_series(1.0, rates)
#' dexp_series(c(0.5, 1.0, 2.0), rates, log = TRUE)
dexp_series <- function(t, rates, log = FALSE) {
    dexp(t, dist.structure::exp_series(rates)$total_rate, log)
}

#' Cumulative distribution function for exponential series.
#'
#' Thin wrapper. Equivalent to
#' `algebraic.dist::cdf(dist.structure::exp_series(rates))(t)`.
#'
#' @param t Vector of series system lifetimes.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param lower.tail Logical; if TRUE (default), probabilities are P(X <= x),
#'                   otherwise, P(X > x).
#' @param log.p Logical; if TRUE, return the log of the cdf.
#' @return The cumulative probabilities evaluated at the specified lifetimes.
#' @note Preserved for backwards compatibility. New code: use
#'   [dist.structure::exp_series()] then [algebraic.dist::cdf()].
#' @importFrom stats pexp
#' @export
#' @examples
#' rates <- c(0.5, 0.3, 0.2)
#' pexp_series(1.0, rates)
#' pexp_series(c(0.5, 1.0, 2.0), rates)
pexp_series <- function(t, rates, lower.tail = TRUE, log.p = FALSE) {
    pexp(t, dist.structure::exp_series(rates)$total_rate, lower.tail, log.p)
}

#' Hazard function for exponential series.
#'
#' Thin wrapper. The series of m independent exponentials has constant
#' system hazard `sum(rates)`. Equivalent to
#' `algebraic.dist::hazard(dist.structure::exp_series(rates))(t)`.
#'
#' @param t Vector of series system lifetimes.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param log.p Logical; if TRUE, return the log of the hazard function.
#' @return The hazard function evaluated at the specified lifetimes.
#' @note Preserved for backwards compatibility. New code: use
#'   [dist.structure::exp_series()] then [algebraic.dist::hazard()].
#' @export
#' @examples
#' rates <- c(0.5, 0.3, 0.2)
#' hazard_exp_series(1.0, rates)
#' hazard_exp_series(c(0.5, 1.0), rates, log.p = TRUE)
hazard_exp_series <- function(t, rates, log.p = FALSE) {
    # The series-of-exponentials hazard is constant = sum(rates) for
    # t >= 0. The original API returned 0 for t < 0 and used a scalar
    # `ifelse(log.p, ...)` wrapper that constrained the output length to
    # 1 (the length of `log.p`); preserved here for backwards
    # compatibility. The dist.structure aggregate is used for the rate
    # value so this remains the canonical source.
    rate_total <- dist.structure::exp_series(rates)$total_rate
    ifelse(log.p,
        ifelse(t < 0, -Inf, log(rate_total)),
        ifelse(t < 0, 0, rate_total)
    )
}


#' Survival function for exponential series.
#'
#' Thin wrapper. Equivalent to
#' `algebraic.dist::surv(dist.structure::exp_series(rates))(t)`.
#'
#' @param t Vector of series system lifetimes.
#' @param rates Vector of rate parameters for exponential component lifetimes.
#' @param log.p Logical; if TRUE, return the log of the survival function.
#' @return The survival function evaluated at the specified lifetimes.
#' @note Preserved for backwards compatibility. New code: use
#'   [dist.structure::exp_series()] then [algebraic.dist::surv()].
#' @importFrom stats pexp
#' @export
#' @examples
#' rates <- c(0.5, 0.3, 0.2)
#' surv.exp_series(1.0, rates)
#' surv.exp_series(c(0.5, 1.0, 2.0), rates)
surv.exp_series <- function(t, rates, log.p = FALSE) {
    pexp(t, dist.structure::exp_series(rates)$total_rate,
         lower.tail = FALSE, log.p = log.p)
}

#' Mean function for exponential series.
#'
#' Computes the expected value of a series system with exponentially
#' distributed component lifetimes. For a series system with component
#' rates \eqn{\lambda_1, \ldots, \lambda_m}, the system lifetime is
#' exponential with rate \eqn{\sum \lambda_j}, so \eqn{E[T] = 1 /
#' \sum \lambda_j}.
#'
#' This method is polymorphic: it accepts both the legacy
#' maskedcauses shape (a numeric vector of rates with class attribute
#' `"exp_series"`) and the current `dist.structure::exp_series()`
#' shape (a list with a `$total_rate` field). This prevents a breaking
#' S3 dispatch collision when both packages are attached.
#'
#' @param x An `exp_series` object. Either a numeric vector of rate
#'   parameters (legacy maskedcauses convention) or a list returned by
#'   [dist.structure::exp_series()].
#' @param ... Additional arguments (ignored, for S3 generic compatibility).
#' @return The mean of the exponential series distribution (1/sum of rates).
#' @note Preserved for backwards compatibility. New code: construct
#'   with `dist.structure::exp_series(rates)` and call `mean()` on it.
#' @method mean exp_series
#' @export
#' @examples
#' # Legacy shape
#' rates <- structure(c(0.5, 0.3, 0.2), class = "exp_series")
#' mean(rates)
#' # Modern shape (dist.structure list) works via the same method
#' mean(dist.structure::exp_series(c(0.5, 0.3, 0.2)))
mean.exp_series <- function(x, ...) {
    # Polymorphic: route both object shapes to the same closed form.
    # The dist.structure list stores the aggregate rate as
    # $total_rate; the legacy numeric vector has the rates directly.
    if (is.list(x) && !is.null(x$total_rate)) {
        return(1 / x$total_rate)
    }
    1 / sum(x)
}
