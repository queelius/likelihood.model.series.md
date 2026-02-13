# Cross-model tests for left-censored and interval-censored observations
#
# Tests cover all three models:
#   - exp_series_md_c1_c2_c3
#   - wei_series_homogeneous_md_c1_c2_c3
#   - wei_series_md_c1_c2_c3
#
# Organized by:
#   1. Hand-computed loglik values (exponential)
#   2. Remark 3: c = {1,...,m} reductions
#   3. Weibull(shape=1) matches exponential for all censoring types
#   4. Score vs numDeriv::grad consistency
#   5. Hessian vs numDeriv::hessian consistency
#   6. Mixed-censoring data
#   7. Edge cases


# =============================================================================
# 1. Hand-computed loglik values — exponential, left-censored
# =============================================================================

test_that("exponential left-censored loglik matches hand computation", {
  # Single left-censored observation: failed before tau=5, C = {1,2}
  # rates = (0.5, 0.3, 0.2), lambda_sys = 1.0, lambda_c = 0.8
  # log L = log(0.8) + log(1 - exp(-1.0 * 5)) - log(1.0)
  #       = log(0.8) + log(1 - exp(-5))

  df <- data.frame(
    t = 5, omega = "left", x1 = TRUE, x2 = TRUE, x3 = FALSE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  rates <- c(0.5, 0.3, 0.2)

  expected <- log(0.8) + log(-expm1(-5))
  expect_equal(ll_fn(df, rates), expected, tolerance = 1e-10)
})

test_that("exponential left-censored loglik, multiple observations", {
  # Two left-censored obs with different taus and candidate sets
  df <- data.frame(
    t = c(3, 8),
    omega = c("left", "left"),
    x1 = c(TRUE, TRUE),
    x2 = c(FALSE, TRUE),
    x3 = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  rates <- c(0.5, 0.3, 0.2)
  lambda_sys <- sum(rates)

  # Obs 1: C={1,3}, lambda_c=0.7
  ll1 <- log(0.7) + log(-expm1(-lambda_sys * 3)) - log(lambda_sys)
  # Obs 2: C={1,2,3}, lambda_c=1.0
  ll2 <- log(1.0) + log(-expm1(-lambda_sys * 8)) - log(lambda_sys)

  expect_equal(ll_fn(df, rates), ll1 + ll2, tolerance = 1e-10)
})


# =============================================================================
# 2. Hand-computed loglik values — exponential, interval-censored
# =============================================================================

test_that("exponential interval-censored loglik matches hand computation", {
  # Single interval obs: failed in (2, 7), C = {1,3}
  # log L = log(lambda_c) - lambda_sys*a + log(1 - exp(-lambda_sys*(b-a))) - log(lambda_sys)
  df <- data.frame(
    t = 2, t_upper = 7, omega = "interval",
    x1 = TRUE, x2 = FALSE, x3 = TRUE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  rates <- c(0.5, 0.3, 0.2)
  lambda_sys <- 1.0
  lambda_c <- 0.7

  expected <- log(lambda_c) - lambda_sys * 2 +
    log(-expm1(-lambda_sys * 5)) - log(lambda_sys)
  expect_equal(ll_fn(df, rates), expected, tolerance = 1e-10)
})


# =============================================================================
# 3. Remark 3: c = {1,...,m} reductions
# =============================================================================

test_that("exp left-censored with c={1,...,m} reduces to log(1-S(tau))", {
  # When all components are in the candidate set, lambda_c = lambda_sys
  # so log L = log(lambda_sys) + log(1-exp(-lambda_sys*tau)) - log(lambda_sys)
  #          = log(1 - exp(-lambda_sys*tau)) = log(1 - S(tau))
  df <- data.frame(
    t = 4, omega = "left", x1 = TRUE, x2 = TRUE, x3 = TRUE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  rates <- c(0.5, 0.3, 0.2)
  lambda_sys <- sum(rates)

  ll_val <- ll_fn(df, rates)
  S_tau <- exp(-lambda_sys * 4)
  expect_equal(ll_val, log(1 - S_tau), tolerance = 1e-10)
})

test_that("exp interval-censored with c={1,...,m} reduces to log(S(a)-S(b))", {
  df <- data.frame(
    t = 2, t_upper = 6, omega = "interval",
    x1 = TRUE, x2 = TRUE, x3 = TRUE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  rates <- c(0.5, 0.3, 0.2)
  lambda_sys <- sum(rates)

  ll_val <- ll_fn(df, rates)
  S_a <- exp(-lambda_sys * 2)
  S_b <- exp(-lambda_sys * 6)
  expect_equal(ll_val, log(S_a - S_b), tolerance = 1e-10)
})

test_that("homogeneous Weibull left-censored with c={1,...,m} reduces to log(1-S(tau))", {
  df <- data.frame(
    t = 50, omega = "left", x1 = TRUE, x2 = TRUE,
    stringsAsFactors = FALSE
  )
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  par <- c(1.5, 100, 200)  # k=1.5, scales=(100,200)

  beta_sys <- wei_series_system_scale(1.5, c(100, 200))
  S_tau <- exp(-(50 / beta_sys)^1.5)
  expected <- log(1 - S_tau)

  expect_equal(ll_fn(df, par), expected, tolerance = 1e-10)
})

test_that("homogeneous Weibull interval-censored with c={1,...,m} reduces to log(S(a)-S(b))", {
  df <- data.frame(
    t = 30, t_upper = 70, omega = "interval",
    x1 = TRUE, x2 = TRUE,
    stringsAsFactors = FALSE
  )
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  par <- c(1.5, 100, 200)

  beta_sys <- wei_series_system_scale(1.5, c(100, 200))
  S_a <- exp(-(30 / beta_sys)^1.5)
  S_b <- exp(-(70 / beta_sys)^1.5)
  expected <- log(S_a - S_b)

  expect_equal(ll_fn(df, par), expected, tolerance = 1e-10)
})


# =============================================================================
# 4. Weibull(shape=1) matches exponential for all censoring types
# =============================================================================

test_that("Weibull(shape=1) left-censored loglik matches exponential", {
  rates <- c(0.5, 0.3, 0.2)
  scales <- 1 / rates

  df <- data.frame(
    t = 5, omega = "left", x1 = TRUE, x2 = TRUE, x3 = FALSE,
    stringsAsFactors = FALSE
  )

  exp_model <- exp_series_md_c1_c2_c3()
  hom_model <- wei_series_homogeneous_md_c1_c2_c3()
  wei_model <- wei_series_md_c1_c2_c3()

  ll_exp <- loglik(exp_model)(df, rates)
  ll_hom <- loglik(hom_model)(df, c(1, scales))
  ll_wei <- loglik(wei_model)(df, as.numeric(rbind(rep(1, 3), scales)))

  expect_equal(ll_hom, ll_exp, tolerance = 1e-8)
  expect_equal(ll_wei, ll_exp, tolerance = 1e-6)
})

test_that("Weibull(shape=1) interval-censored loglik matches exponential", {
  rates <- c(0.5, 0.3, 0.2)
  scales <- 1 / rates

  df <- data.frame(
    t = 2, t_upper = 7, omega = "interval",
    x1 = TRUE, x2 = FALSE, x3 = TRUE,
    stringsAsFactors = FALSE
  )

  exp_model <- exp_series_md_c1_c2_c3()
  hom_model <- wei_series_homogeneous_md_c1_c2_c3()
  wei_model <- wei_series_md_c1_c2_c3()

  ll_exp <- loglik(exp_model)(df, rates)
  ll_hom <- loglik(hom_model)(df, c(1, scales))
  ll_wei <- loglik(wei_model)(df, as.numeric(rbind(rep(1, 3), scales)))

  expect_equal(ll_hom, ll_exp, tolerance = 1e-8)
  expect_equal(ll_wei, ll_exp, tolerance = 1e-6)
})


# =============================================================================
# 5. Score vs numDeriv::grad consistency — all models, all censoring types
# =============================================================================

test_that("exponential score matches numDeriv::grad for left-censored", {
  df <- data.frame(
    t = c(3, 8),
    omega = c("left", "left"),
    x1 = c(TRUE, TRUE), x2 = c(FALSE, TRUE), x3 = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3, 0.2)

  analytical <- score(model)(df, rates)
  numerical <- numDeriv::grad(loglik(model), x = rates, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-6)
})

test_that("exponential score matches numDeriv::grad for interval-censored", {
  df <- data.frame(
    t = c(2, 5), t_upper = c(7, 12),
    omega = c("interval", "interval"),
    x1 = c(TRUE, TRUE), x2 = c(FALSE, TRUE), x3 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3, 0.2)

  analytical <- score(model)(df, rates)
  numerical <- numDeriv::grad(loglik(model), x = rates, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-6)
})

test_that("exponential score matches numDeriv::grad for mixed censoring", {
  df <- data.frame(
    t = c(4, 10, 3, 2),
    t_upper = c(NA, NA, NA, 8),
    omega = c("exact", "right", "left", "interval"),
    x1 = c(TRUE, FALSE, TRUE, TRUE),
    x2 = c(TRUE, FALSE, FALSE, TRUE),
    x3 = c(FALSE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3, 0.2)

  analytical <- score(model)(df, rates)
  numerical <- numDeriv::grad(loglik(model), x = rates, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-6)
})

test_that("homogeneous Weibull score matches numDeriv::grad for left-censored", {
  df <- data.frame(
    t = c(30, 60),
    omega = c("left", "left"),
    x1 = c(TRUE, TRUE), x2 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- wei_series_homogeneous_md_c1_c2_c3()
  par <- c(1.5, 100, 200)

  analytical <- score(model)(df, par)
  numerical <- numDeriv::grad(loglik(model), x = par, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-5)
})

test_that("homogeneous Weibull score matches numDeriv::grad for interval-censored", {
  df <- data.frame(
    t = c(30, 40), t_upper = c(60, 80),
    omega = c("interval", "interval"),
    x1 = c(TRUE, TRUE), x2 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- wei_series_homogeneous_md_c1_c2_c3()
  par <- c(1.5, 100, 200)

  analytical <- score(model)(df, par)
  numerical <- numDeriv::grad(loglik(model), x = par, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-5)
})

test_that("heterogeneous Weibull score matches numDeriv::grad for left-censored", {
  df <- data.frame(
    t = c(30, 60),
    omega = c("left", "left"),
    x1 = c(TRUE, TRUE), x2 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- wei_series_md_c1_c2_c3()
  par <- c(1.5, 100, 2.0, 200)

  analytical <- score(model)(df, par)
  numerical <- numDeriv::grad(loglik(model), x = par, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-5)
})

test_that("heterogeneous Weibull score matches numDeriv::grad for interval-censored", {
  df <- data.frame(
    t = c(30, 40), t_upper = c(60, 80),
    omega = c("interval", "interval"),
    x1 = c(TRUE, TRUE), x2 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- wei_series_md_c1_c2_c3()
  par <- c(1.5, 100, 2.0, 200)

  analytical <- score(model)(df, par)
  numerical <- numDeriv::grad(loglik(model), x = par, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-5)
})


# =============================================================================
# 6. Hessian vs numDeriv::hessian consistency
# =============================================================================

test_that("exponential hessian matches numDeriv::hessian for left-censored", {
  df <- data.frame(
    t = c(3, 8),
    omega = c("left", "left"),
    x1 = c(TRUE, TRUE), x2 = c(FALSE, TRUE), x3 = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3, 0.2)

  analytical <- hess_loglik(model)(df, rates)
  numerical <- numDeriv::hessian(loglik(model), x = rates, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-5)
})

test_that("exponential hessian matches numDeriv::hessian for interval-censored", {
  df <- data.frame(
    t = c(2, 5), t_upper = c(7, 12),
    omega = c("interval", "interval"),
    x1 = c(TRUE, TRUE), x2 = c(FALSE, TRUE), x3 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3, 0.2)

  analytical <- hess_loglik(model)(df, rates)
  numerical <- numDeriv::hessian(loglik(model), x = rates, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-5)
})

test_that("exponential hessian matches numDeriv::hessian for mixed censoring", {
  df <- data.frame(
    t = c(4, 10, 3, 2),
    t_upper = c(NA, NA, NA, 8),
    omega = c("exact", "right", "left", "interval"),
    x1 = c(TRUE, FALSE, TRUE, TRUE),
    x2 = c(TRUE, FALSE, FALSE, TRUE),
    x3 = c(FALSE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3, 0.2)

  analytical <- hess_loglik(model)(df, rates)
  numerical <- numDeriv::hessian(loglik(model), x = rates, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-5)
})

test_that("homogeneous Weibull hessian matches numDeriv::hessian for left+interval", {
  df <- data.frame(
    t = c(30, 30), t_upper = c(NA, 70),
    omega = c("left", "interval"),
    x1 = c(TRUE, TRUE), x2 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- wei_series_homogeneous_md_c1_c2_c3()
  par <- c(1.5, 100, 200)

  analytical <- hess_loglik(model)(df, par)
  numerical <- numDeriv::hessian(loglik(model), x = par, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-4)
})

test_that("heterogeneous Weibull hessian matches numDeriv::hessian for left+interval", {
  df <- data.frame(
    t = c(30, 30), t_upper = c(NA, 70),
    omega = c("left", "interval"),
    x1 = c(TRUE, TRUE), x2 = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  model <- wei_series_md_c1_c2_c3()
  par <- c(1.5, 100, 2.0, 200)

  analytical <- hess_loglik(model)(df, par)
  numerical <- numDeriv::hessian(loglik(model), x = par, df = df)

  expect_equal(analytical, numerical, tolerance = 1e-4)
})


# =============================================================================
# 7. Mixed-censoring data — all four types together
# =============================================================================

test_that("loglik is finite for mixed-censoring data across all models", {
  df <- data.frame(
    t = c(5, 15, 4, 3),
    t_upper = c(NA, NA, NA, 10),
    omega = c("exact", "right", "left", "interval"),
    x1 = c(TRUE, FALSE, TRUE, TRUE),
    x2 = c(TRUE, FALSE, TRUE, FALSE),
    x3 = c(FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  # Exponential
  exp_model <- exp_series_md_c1_c2_c3()
  ll_exp <- loglik(exp_model)(df, c(0.5, 0.3, 0.2))
  expect_true(is.finite(ll_exp))

  # Homogeneous Weibull
  hom_model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_hom <- loglik(hom_model)(df, c(1.5, 100, 200, 300))
  expect_true(is.finite(ll_hom))

  # Heterogeneous Weibull
  wei_model <- wei_series_md_c1_c2_c3()
  ll_wei <- loglik(wei_model)(df, c(1.5, 100, 2.0, 200, 1.2, 300))
  expect_true(is.finite(ll_wei))
})

test_that("score returns correct-length vector for mixed-censoring data", {
  df <- data.frame(
    t = c(5, 15, 4, 3),
    t_upper = c(NA, NA, NA, 10),
    omega = c("exact", "right", "left", "interval"),
    x1 = c(TRUE, FALSE, TRUE, TRUE),
    x2 = c(TRUE, FALSE, TRUE, FALSE),
    x3 = c(FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  exp_model <- exp_series_md_c1_c2_c3()
  s_exp <- score(exp_model)(df, c(0.5, 0.3, 0.2))
  expect_length(s_exp, 3)
  expect_true(all(is.finite(s_exp)))

  hom_model <- wei_series_homogeneous_md_c1_c2_c3()
  s_hom <- score(hom_model)(df, c(1.5, 100, 200, 300))
  expect_length(s_hom, 4)
  expect_true(all(is.finite(s_hom)))

  wei_model <- wei_series_md_c1_c2_c3()
  s_wei <- score(wei_model)(df, c(1.5, 100, 2.0, 200, 1.2, 300))
  expect_length(s_wei, 6)
  expect_true(all(is.finite(s_wei)))
})


# =============================================================================
# 8. Edge cases
# =============================================================================

test_that("all observations left-censored works", {
  df <- data.frame(
    t = c(5, 8, 3),
    omega = c("left", "left", "left"),
    x1 = c(TRUE, TRUE, TRUE),
    x2 = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  model <- exp_series_md_c1_c2_c3()
  ll <- loglik(model)(df, c(0.5, 0.3))
  s <- score(model)(df, c(0.5, 0.3))
  h <- hess_loglik(model)(df, c(0.5, 0.3))

  expect_true(is.finite(ll))
  expect_length(s, 2)
  expect_true(all(is.finite(s)))
  expect_equal(dim(h), c(2, 2))
  expect_true(all(is.finite(h)))
})

test_that("all observations interval-censored works", {
  df <- data.frame(
    t = c(2, 5, 1),
    t_upper = c(6, 10, 4),
    omega = c("interval", "interval", "interval"),
    x1 = c(TRUE, TRUE, TRUE),
    x2 = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  model <- exp_series_md_c1_c2_c3()
  ll <- loglik(model)(df, c(0.5, 0.3))
  s <- score(model)(df, c(0.5, 0.3))
  h <- hess_loglik(model)(df, c(0.5, 0.3))

  expect_true(is.finite(ll))
  expect_length(s, 2)
  expect_true(all(is.finite(s)))
  expect_equal(dim(h), c(2, 2))
  expect_true(all(is.finite(h)))
})

test_that("single observation of each type works", {
  # Exact
  df_exact <- data.frame(
    t = 5, omega = "exact", x1 = TRUE, x2 = FALSE,
    stringsAsFactors = FALSE
  )
  # Right
  df_right <- data.frame(
    t = 10, omega = "right", x1 = FALSE, x2 = FALSE,
    stringsAsFactors = FALSE
  )
  # Left
  df_left <- data.frame(
    t = 3, omega = "left", x1 = TRUE, x2 = TRUE,
    stringsAsFactors = FALSE
  )
  # Interval
  df_int <- data.frame(
    t = 2, t_upper = 8, omega = "interval", x1 = TRUE, x2 = FALSE,
    stringsAsFactors = FALSE
  )

  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3)

  for (df_obs in list(df_exact, df_right, df_left, df_int)) {
    ll <- loglik(model)(df_obs, rates)
    expect_true(is.finite(ll))
    s <- score(model)(df_obs, rates)
    expect_length(s, 2)
    expect_true(all(is.finite(s)))
  }
})

test_that("interval-censored validation: t_upper must exceed t", {
  df <- data.frame(
    t = 5, t_upper = 3, omega = "interval", x1 = TRUE, x2 = FALSE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  expect_error(loglik(model)(df, c(0.5, 0.3)), "t_upper > t")
})

test_that("left-censored must have non-empty candidate set", {
  df <- data.frame(
    t = 5, omega = "left", x1 = FALSE, x2 = FALSE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  expect_error(loglik(model)(df, c(0.5, 0.3)), "non-empty candidate set")
})

test_that("interval-censored must have non-empty candidate set", {
  df <- data.frame(
    t = 2, t_upper = 5, omega = "interval", x1 = FALSE, x2 = FALSE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  expect_error(loglik(model)(df, c(0.5, 0.3)), "non-empty candidate set")
})

test_that("interval-censored requires t_upper column", {
  df <- data.frame(
    t = 2, omega = "interval", x1 = TRUE, x2 = FALSE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  expect_error(loglik(model)(df, c(0.5, 0.3)), "t_upper")
})

test_that("invalid omega values are rejected", {
  df <- data.frame(
    t = 5, omega = "unknown", x1 = TRUE, x2 = FALSE,
    stringsAsFactors = FALSE
  )
  model <- exp_series_md_c1_c2_c3()
  expect_error(loglik(model)(df, c(0.5, 0.3)), "invalid omega")
})


# =============================================================================
# 9. Heterogeneous Weibull integrand correctness
# =============================================================================

test_that("het Weibull left-censored with c={1,...,m} reduces to 1-S(tau)", {
  # When c contains all components, the integral of sum_j h_j(t)*S(t) from 0 to tau
  # equals the CDF F(tau) = 1 - S(tau)
  df <- data.frame(
    t = 50, omega = "left", x1 = TRUE, x2 = TRUE,
    stringsAsFactors = FALSE
  )
  model <- wei_series_md_c1_c2_c3()
  par <- c(1.5, 100, 2.0, 200)

  ll_val <- loglik(model)(df, par)
  # S(tau) = exp(-sum_j (tau/beta_j)^k_j)
  S_tau <- exp(-sum((50 / c(100, 200))^c(1.5, 2.0)))
  expected <- log(1 - S_tau)

  expect_equal(ll_val, expected, tolerance = 1e-6)
})

test_that("het Weibull interval-censored with c={1,...,m} reduces to log(S(a)-S(b))", {
  df <- data.frame(
    t = 30, t_upper = 70, omega = "interval",
    x1 = TRUE, x2 = TRUE,
    stringsAsFactors = FALSE
  )
  model <- wei_series_md_c1_c2_c3()
  par <- c(1.5, 100, 2.0, 200)

  ll_val <- loglik(model)(df, par)
  S_a <- exp(-sum((30 / c(100, 200))^c(1.5, 2.0)))
  S_b <- exp(-sum((70 / c(100, 200))^c(1.5, 2.0)))
  expected <- log(S_a - S_b)

  expect_equal(ll_val, expected, tolerance = 1e-6)
})
