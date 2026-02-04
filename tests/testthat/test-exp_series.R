# Tests for exponential series distribution functions

# ==============================================================================
# Tests for rexp_series
# ==============================================================================

test_that("rexp_series returns vector of correct length", {
  n <- 100
  rates <- c(0.5, 0.3, 0.2)

  x <- rexp_series(n, rates)

  expect_length(x, n)
  expect_true(all(x >= 0))
  expect_type(x, "double")
})

test_that("rexp_series with keep_latent returns matrix with correct dimensions", {
  n <- 100
  rates <- c(0.5, 0.3, 0.2)
  m <- length(rates)

  x <- rexp_series(n, rates, keep_latent = TRUE)

  expect_true(is.matrix(x))
  expect_equal(dim(x), c(n, m + 1))
  expect_true(all(x >= 0))
})

test_that("rexp_series with keep_latent has system lifetime as minimum", {
  set.seed(42)
  n <- 50
  rates <- c(0.5, 0.3, 0.2)

  x <- rexp_series(n, rates, keep_latent = TRUE)

  # First column is system lifetime (minimum of component lifetimes)
  system_lifetime <- x[, 1]
  component_lifetimes <- x[, -1]
  min_comp <- apply(component_lifetimes, 1, min)

  expect_equal(system_lifetime, min_comp)
})

test_that("rexp_series sample mean converges to theoretical mean", {
  set.seed(123)
  n <- 10000
  rates <- c(0.5, 0.3, 0.2)

  x <- rexp_series(n, rates)
  sample_mean <- mean(x)

  # For exponential series, rate = sum(rates), mean = 1/rate
  theoretical_mean <- 1 / sum(rates)

  expect_equal(sample_mean, theoretical_mean, tolerance = 0.05)
})

test_that("rexp_series handles single component", {
  n <- 100
  rates <- 0.5

  x <- rexp_series(n, rates)

  expect_length(x, n)
  expect_true(all(x >= 0))
})


# ==============================================================================
# Tests for dexp_series
# ==============================================================================

test_that("dexp_series returns correct density values", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(0, 1, 2, 5)

  d <- dexp_series(t, rates)

  # Exponential series has rate = sum(rates) = 1.0
  expected <- dexp(t, rate = sum(rates))

  expect_equal(d, expected)
})

test_that("dexp_series returns log density when log = TRUE", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(0.5, 1, 2)

  log_d <- dexp_series(t, rates, log = TRUE)
  d <- dexp_series(t, rates, log = FALSE)

  expect_equal(log_d, log(d))
})

test_that("dexp_series integrates to 1", {
  rates <- c(0.5, 0.3, 0.2)

  # Numerical integration
  result <- integrate(function(x) dexp_series(x, rates), 0, Inf)

  expect_equal(result$value, 1, tolerance = 1e-6)
})

test_that("dexp_series returns 0 for negative values", {
  rates <- c(0.5, 0.3, 0.2)

  d <- dexp_series(-1, rates)

  expect_equal(d, 0)
})


# ==============================================================================
# Tests for pexp_series
# ==============================================================================

test_that("pexp_series returns correct CDF values", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(0, 1, 2, 5)

  p <- pexp_series(t, rates)

  # Exponential series has rate = sum(rates) = 1.0
  expected <- pexp(t, rate = sum(rates))

  expect_equal(p, expected)
})

test_that("pexp_series returns correct values with lower.tail = FALSE", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(1, 2, 5)

  p_lower <- pexp_series(t, rates, lower.tail = TRUE)
  p_upper <- pexp_series(t, rates, lower.tail = FALSE)

  expect_equal(p_lower + p_upper, rep(1, length(t)))
})

test_that("pexp_series returns log probabilities when log.p = TRUE", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(1, 2, 5)

  log_p <- pexp_series(t, rates, log.p = TRUE)
  p <- pexp_series(t, rates, log.p = FALSE)

  expect_equal(log_p, log(p))
})

test_that("pexp_series has correct boundary behavior", {
  rates <- c(0.5, 0.3, 0.2)

  expect_equal(pexp_series(0, rates), 0)
  expect_equal(pexp_series(Inf, rates), 1)
})

test_that("pexp_series is consistent with dexp_series", {
  rates <- c(0.5, 0.3, 0.2)
  t_vals <- c(0.5, 1, 2, 5)

  for (t in t_vals) {
    # CDF should be integral of PDF from 0 to t
    result <- integrate(function(x) dexp_series(x, rates), 0, t)
    cdf_val <- pexp_series(t, rates)

    expect_equal(result$value, cdf_val, tolerance = 1e-6)
  }
})


# ==============================================================================
# Tests for qexp_series
# ==============================================================================

test_that("qexp_series returns correct quantile values", {
  rates <- c(0.5, 0.3, 0.2)
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)

  q <- qexp_series(p, rates)

  # Exponential series has rate = sum(rates)
  expected <- qexp(p, rate = sum(rates))

  expect_equal(q, expected)
})

test_that("qexp_series is inverse of pexp_series", {
  rates <- c(0.5, 0.3, 0.2)
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)

  q <- qexp_series(p, rates)
  p_back <- pexp_series(q, rates)

  expect_equal(p, p_back, tolerance = 1e-10)
})

test_that("qexp_series handles boundary probabilities", {
  rates <- c(0.5, 0.3, 0.2)

  expect_equal(qexp_series(0, rates), 0)
  expect_equal(qexp_series(1, rates), Inf)
})

test_that("qexp_series with lower.tail = FALSE", {
  rates <- c(0.5, 0.3, 0.2)
  p <- 0.25

  q_lower <- qexp_series(p, rates, lower.tail = TRUE)
  q_upper <- qexp_series(1 - p, rates, lower.tail = FALSE)

  expect_equal(q_lower, q_upper, tolerance = 1e-10)
})

test_that("qexp_series with log.p = TRUE", {
  rates <- c(0.5, 0.3, 0.2)
  p <- 0.25

  q <- qexp_series(p, rates, log.p = FALSE)
  q_log <- qexp_series(log(p), rates, log.p = TRUE)

  expect_equal(q, q_log, tolerance = 1e-10)
})


# ==============================================================================
# Tests for hazard_exp_series
# ==============================================================================

test_that("hazard_exp_series returns constant hazard rate", {
  rates <- c(0.5, 0.3, 0.2)
  t <- 1

  h <- hazard_exp_series(t, rates)

  # Hazard for exponential is constant = rate = sum(rates)
  # Note: current implementation returns scalar, not vector
  expect_equal(h, sum(rates))
})

test_that("hazard_exp_series with log.p = TRUE returns log of hazard rate", {
  rates <- c(0.5, 0.3, 0.2)
  t <- 1

  h <- hazard_exp_series(t, rates, log.p = FALSE)
  log_h <- hazard_exp_series(t, rates, log.p = TRUE)

  expect_equal(h, sum(rates))
  expect_equal(log_h, log(sum(rates)))
})

test_that("hazard_exp_series returns 0 for negative times", {
  rates <- c(0.5, 0.3, 0.2)

  h <- hazard_exp_series(-1, rates, log.p = FALSE)

  expect_equal(h, 0)
})

test_that("hazard_exp_series returns -Inf for negative times with log.p = TRUE", {
  rates <- c(0.5, 0.3, 0.2)

  h <- hazard_exp_series(-1, rates, log.p = TRUE)

  expect_equal(h, -Inf)
})


# ==============================================================================
# Tests for surv.exp_series
# ==============================================================================

test_that("surv.exp_series returns correct survival function S(t) = P(T > t)", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(0, 1, 2, 5)

  s <- surv.exp_series(t, rates)

  # Survival function S(t) = 1 - F(t) = exp(-lambda * t)
  expected <- pexp(t, sum(rates), lower.tail = FALSE)

  expect_equal(s, expected)
})

test_that("surv.exp_series returns log survival when log.p = TRUE", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(1, 2, 5)

  s <- surv.exp_series(t, rates, log.p = FALSE)
  log_s <- surv.exp_series(t, rates, log.p = TRUE)

  expect_equal(log_s, log(s), tolerance = 1e-10)
})

test_that("surv.exp_series has correct boundary behavior", {
  rates <- c(0.5, 0.3, 0.2)

  # S(0) = 1, S(Inf) = 0
  expect_equal(surv.exp_series(0, rates), 1)
  expect_equal(surv.exp_series(Inf, rates), 0)
})

test_that("surv.exp_series + pexp_series = 1 (survival + CDF)", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(0.5, 1, 2, 5)

  s <- surv.exp_series(t, rates)
  p <- pexp_series(t, rates)

  expect_equal(s + p, rep(1, length(t)), tolerance = 1e-10)
})


# ==============================================================================
# Tests for mean.exp_series
# ==============================================================================

test_that("mean.exp_series returns correct mean", {
  rates <- c(0.5, 0.3, 0.2)
  class(rates) <- "exp_series"

  m <- mean(rates)

  # Mean of exponential with rate lambda is 1/lambda
  expected <- 1 / sum(rates)

  expect_equal(m, expected)
})

test_that("mean.exp_series handles single component", {
  rates <- 0.5
  class(rates) <- "exp_series"

  m <- mean(rates)

  expect_equal(m, 1 / 0.5)
})

test_that("mean.exp_series is consistent with sample mean from rexp_series", {
  set.seed(42)
  rates <- c(0.5, 0.3, 0.2)
  class(rates) <- "exp_series"

  theoretical_mean <- mean(rates)
  sample <- rexp_series(50000, rates)
  sample_mean <- mean(sample)

  expect_equal(sample_mean, theoretical_mean, tolerance = 0.02)
})


# ==============================================================================
# Edge case tests
# ==============================================================================

test_that("distribution functions handle equal rates", {
  rates <- c(0.5, 0.5, 0.5)
  t <- 1

  d <- dexp_series(t, rates)
  p <- pexp_series(t, rates)
  q <- qexp_series(0.5, rates)
  h <- hazard_exp_series(t, rates)
  s <- surv.exp_series(t, rates)

  expect_true(all(is.finite(c(d, p, q, h, s))))
})

test_that("distribution functions handle very small rates", {
  rates <- c(0.001, 0.002, 0.003)
  t <- 100

  d <- dexp_series(t, rates)
  p <- pexp_series(t, rates)
  s <- surv.exp_series(t, rates)

  expect_true(all(is.finite(c(d, p, s))))
  expect_true(d > 0)
  expect_true(p > 0 && p < 1)
  expect_true(s > 0 && s < 1)
})

test_that("distribution functions handle very large rates", {
  rates <- c(10, 20, 30)
  t <- 0.01

  d <- dexp_series(t, rates)
  p <- pexp_series(t, rates)
  s <- surv.exp_series(t, rates)

  expect_true(all(is.finite(c(d, p, s))))
})

test_that("vectorized operations work correctly for d, p, s functions", {
  rates <- c(0.5, 0.3, 0.2)
  t <- c(0, 0.5, 1, 2, 5, 10)

  d <- dexp_series(t, rates)
  p <- pexp_series(t, rates)
  s <- surv.exp_series(t, rates)

  expect_length(d, 6)
  expect_length(p, 6)
  expect_length(s, 6)

  # All should be non-negative
  expect_true(all(d >= 0))
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(s >= 0 & s <= 1))
})

test_that("hazard_exp_series returns scalar constant", {
  # Note: hazard_exp_series returns a single scalar (the constant hazard rate)
  # rather than vectorizing over t
  rates <- c(0.5, 0.3, 0.2)
  t <- c(0, 1, 2)

  h <- hazard_exp_series(t, rates)

  # Should return the constant hazard = sum(rates)
  expect_equal(h, sum(rates))
})
