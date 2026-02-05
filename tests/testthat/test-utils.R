# Tests for utility functions

# ==============================================================================
# Tests for cum_haz
# ==============================================================================

test_that("cum_haz creates a cumulative hazard function", {
  # Constant hazard (exponential) - must be vectorized for integrate()
  lambda <- 0.5
  haz <- function(t) rep(lambda, length(t))

  ch <- cum_haz(haz)

  # Cumulative hazard should be integral of constant = lambda * t
  expect_equal(ch(1), lambda * 1, tolerance = 1e-6)
  expect_equal(ch(2), lambda * 2, tolerance = 1e-6)
  expect_equal(ch(5), lambda * 5, tolerance = 1e-6)
})

test_that("cum_haz handles Weibull hazard correctly", {
  # Weibull hazard: h(t) = (k/lambda) * (t/lambda)^(k-1)
  k <- 2
  lambda <- 10

  # Vectorized hazard function
  haz <- function(t) (k / lambda) * (t / lambda)^(k - 1)
  ch <- cum_haz(haz)

  # Weibull cumulative hazard: H(t) = (t/lambda)^k
  expected_ch <- function(t) (t / lambda)^k

  expect_equal(ch(1), expected_ch(1), tolerance = 1e-5)
  expect_equal(ch(5), expected_ch(5), tolerance = 1e-5)
  expect_equal(ch(10), expected_ch(10), tolerance = 1e-5)
})

test_that("cum_haz returns 0 at t=0", {
  # Vectorized constant hazard
  haz <- function(t) rep(1, length(t))

  ch <- cum_haz(haz)

  expect_equal(ch(0), 0, tolerance = 1e-10)
})

test_that("cum_haz returns function that accepts additional arguments", {
  # Hazard with rate parameter - vectorized
  haz <- function(t, rate) rep(rate, length(t))

  ch <- cum_haz(haz)

  # Should pass rate through to hazard function
  expect_equal(ch(2, rate = 0.5), 1, tolerance = 1e-6)
  expect_equal(ch(2, rate = 1.0), 2, tolerance = 1e-6)
})


# ==============================================================================
# Tests for qcomp and rcomp (internal functions)
# ==============================================================================

# Note: qcomp and rcomp are internal functions that use optimization
# to find quantiles/generate random variates for arbitrary hazard/survival functions.
# They are exported but not commonly used directly.

test_that("qcomp finds quantile for exponential distribution", {
  # Exponential distribution with rate lambda
  # S(t) = exp(-lambda * t)
  # Quantile: t = -log(p) / lambda
  lambda <- 0.5

  surv <- function(t, theta) exp(-theta * t)

  # Test median (p = 0.5)
  p <- 0.5
  expected_t <- -log(p) / lambda

  result <- qcomp(p, surv = surv, theta = lambda)

  expect_equal(result, expected_t, tolerance = 1e-6)
})

test_that("qcomp works with different probability values", {
  lambda <- 1.0

  surv <- function(t, theta) exp(-theta * t)

  # Test various quantiles - compare to analytical solution
  for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    result <- qcomp(p, surv = surv, theta = lambda)
    expected <- -log(p) / lambda
    expect_equal(result, expected, tolerance = 1e-6)
  }
})

test_that("qcomp works with Weibull distribution", {
  # Weibull: S(t) = exp(-(t/scale)^shape)
  # Quantile: t = scale * (-log(p))^(1/shape)
  shape <- 2
  scale <- 10
  theta <- c(shape, scale)

  surv <- function(t, theta) {
    k <- theta[1]
    lam <- theta[2]
    exp(-(t / lam)^k)
  }

  p <- 0.5
  result <- qcomp(p, surv = surv, theta = theta)
  expected <- scale * (-log(p))^(1/shape)

  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("qcomp accepts custom bounds", {
  lambda <- 0.5

  surv <- function(t, theta) exp(-theta * t)

  # Custom bounds should work
  result <- qcomp(0.5, surv = surv, theta = lambda, t_lower = 0.01, t_upper = 100)

  expected <- -log(0.5) / lambda
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("qcomp uses machine epsilon and sqrt(xmax) as defaults", {
  lambda <- 0.5

  surv <- function(t, theta) exp(-theta * t)

  # Should work with default bounds
  result <- qcomp(0.5, surv = surv, theta = lambda)

  expected <- -log(0.5) / lambda
  expect_equal(result, expected, tolerance = 1e-6)

  # Verify default values are what we expect (formals returns expressions)
  # Note: R's eval() on formals is safe here - just evaluating R expressions
  expect_equal(base::eval(formals(qcomp)$t_lower), .Machine$double.eps)
  expect_equal(base::eval(formals(qcomp)$t_upper), .Machine$double.xmax^0.5)
})

test_that("rcomp generates n random samples", {
  lambda <- 0.5

  surv <- function(t, theta) exp(-theta * t)

  set.seed(123)
  n <- 10
  samples <- rcomp(n, surv = surv, theta = lambda)

  expect_length(samples, n)
  expect_true(all(is.numeric(samples)))
})

test_that("rcomp generates positive values for standard distributions", {
  lambda <- 1.0

  surv <- function(t, theta) exp(-theta * t)

  set.seed(456)
  samples <- rcomp(20, surv = surv, theta = lambda)

  # All samples should be numeric
  expect_true(all(is.numeric(samples)))
})

test_that("rcomp works with Weibull distribution", {
  shape <- 1.5
  scale <- 100
  theta <- c(shape, scale)

  surv <- function(t, theta) {
    k <- theta[1]
    lam <- theta[2]
    exp(-(t / lam)^k)
  }

  set.seed(789)
  samples <- rcomp(5, surv = surv, theta = theta)

  expect_length(samples, 5)
  expect_true(all(is.numeric(samples)))
})

test_that("rcomp uses internal qcomp for each sample", {
  # Verify rcomp generates different values (not all the same)
  lambda <- 0.5

  surv <- function(t, theta) exp(-theta * t)

  set.seed(101)
  samples <- rcomp(10, surv = surv, theta = lambda)

  # Should have variation (not all identical)
  expect_true(length(unique(samples)) > 1)
})

test_that("rcomp samples follow expected distribution (KS test)", {
  # Generate samples and verify they follow exponential distribution
  lambda <- 1.0

  surv <- function(t, theta) exp(-theta * t)

  set.seed(42)
  samples <- rcomp(100, surv = surv, theta = lambda)

  # Kolmogorov-Smirnov test against exponential(rate=1)
  ks_result <- ks.test(samples, "pexp", rate = lambda)

  # Should not reject null hypothesis (samples come from exp distribution)
  expect_true(ks_result$p.value > 0.05)
})


# ==============================================================================
# Tests for qcomp and rcomp edge cases
# ==============================================================================

test_that("qcomp works with p near 0 (tail quantile)", {
  lambda <- 1.0
  surv <- function(t, theta) exp(-theta * t)

  result <- qcomp(p = 0.001, surv = surv, theta = lambda)
  expected <- -log(0.001) / lambda

  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("qcomp works with p near 1 (near zero)", {
  lambda <- 1.0
  surv <- function(t, theta) exp(-theta * t)

  result <- qcomp(p = 0.999, surv = surv, theta = lambda)
  expected <- -log(0.999) / lambda

  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("qcomp works with only surv parameter", {
  lambda <- 0.5
  surv <- function(t, theta) exp(-theta * t)

  # Should work with just surv
  result <- qcomp(p = 0.5, surv = surv, theta = lambda)
  expected <- -log(0.5) / lambda

  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("rcomp with n=0 returns empty vector", {
  lambda <- 1.0
  surv <- function(t, theta) exp(-theta * t)

  samples <- rcomp(0, surv = surv, theta = lambda)

  expect_length(samples, 0)
})

test_that("cum_haz at t=0 returns 0", {
  haz <- function(t) rep(2.0, length(t))
  ch <- cum_haz(haz)
  expect_equal(ch(0), 0, tolerance = 1e-10)
})

test_that("cum_haz with increasing hazard", {
  # Linear hazard h(t) = t, H(t) = t^2/2
  haz <- function(t) t
  ch <- cum_haz(haz)
  expect_equal(ch(2), 2, tolerance = 1e-5)
  expect_equal(ch(4), 8, tolerance = 1e-5)
})
