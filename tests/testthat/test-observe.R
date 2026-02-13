# Tests for observation functor system (R/observe.R)

# =============================================================================
# 1. observe_right_censor
# =============================================================================

test_that("observe_right_censor returns exact for t_true <= tau", {
  obs <- observe_right_censor(tau = 10)
  result <- obs(5)

  expect_equal(result$t, 5)
  expect_equal(result$omega, "exact")
  expect_true(is.na(result$t_upper))
})

test_that("observe_right_censor returns right for t_true > tau", {
  obs <- observe_right_censor(tau = 10)
  result <- obs(15)

  expect_equal(result$t, 10)
  expect_equal(result$omega, "right")
  expect_true(is.na(result$t_upper))
})

test_that("observe_right_censor boundary: t_true == tau is exact", {
  obs <- observe_right_censor(tau = 10)
  result <- obs(10)

  expect_equal(result$t, 10)
  expect_equal(result$omega, "exact")
})

test_that("observe_right_censor with tau=Inf: all exact", {
  obs <- observe_right_censor(tau = Inf)
  result <- obs(1e6)

  expect_equal(result$t, 1e6)
  expect_equal(result$omega, "exact")
})


# =============================================================================
# 2. observe_left_censor
# =============================================================================

test_that("observe_left_censor returns left for t_true <= tau", {
  obs <- observe_left_censor(tau = 10)
  result <- obs(5)

  expect_equal(result$t, 10)
  expect_equal(result$omega, "left")
  expect_true(is.na(result$t_upper))
})

test_that("observe_left_censor returns right for t_true > tau", {
  obs <- observe_left_censor(tau = 10)
  result <- obs(15)

  expect_equal(result$t, 10)
  expect_equal(result$omega, "right")
  expect_true(is.na(result$t_upper))
})

test_that("observe_left_censor boundary: t_true == tau is left", {
  obs <- observe_left_censor(tau = 10)
  result <- obs(10)

  expect_equal(result$t, 10)
  expect_equal(result$omega, "left")
})


# =============================================================================
# 3. observe_periodic
# =============================================================================

test_that("observe_periodic produces interval observations", {
  obs <- observe_periodic(delta = 10, tau = 100)
  result <- obs(25)

  expect_equal(result$t, 20)
  expect_equal(result$omega, "interval")
  expect_equal(result$t_upper, 30)
})

test_that("observe_periodic: failure at exact inspection boundary", {
  obs <- observe_periodic(delta = 10, tau = 100)
  result <- obs(30)

  expect_equal(result$t, 30)
  expect_equal(result$omega, "interval")
  expect_equal(result$t_upper, 40)
})

test_that("observe_periodic: failure in first interval", {
  obs <- observe_periodic(delta = 5, tau = 100)
  result <- obs(3)

  expect_equal(result$t, 0)
  expect_equal(result$omega, "interval")
  expect_equal(result$t_upper, 5)
})

test_that("observe_periodic: right-censored when t_true > tau", {
  obs <- observe_periodic(delta = 10, tau = 100)
  result <- obs(150)

  expect_equal(result$t, 100)
  expect_equal(result$omega, "right")
  expect_true(is.na(result$t_upper))
})

test_that("observe_periodic: tau=Inf means no right-censoring", {
  obs <- observe_periodic(delta = 10, tau = Inf)
  result <- obs(9999)

  expect_equal(result$omega, "interval")
  expect_equal(result$t, 9990)
  expect_equal(result$t_upper, 10000)
})


# =============================================================================
# 4. observe_mixture
# =============================================================================

test_that("observe_mixture produces valid observations", {
  set.seed(42)
  obs <- observe_mixture(
    observe_right_censor(tau = 100),
    observe_left_censor(tau = 50),
    weights = c(0.5, 0.5)
  )

  # Generate many observations and check all are valid
  omegas <- character(100)
  for (i in seq_len(100)) {
    result <- obs(30)
    omegas[i] <- result$omega
  }

  # Should see both "exact" (from right_censor, since 30 < 100) and "left"
  # (from left_censor, since 30 < 50)
  expect_true("exact" %in% omegas)
  expect_true("left" %in% omegas)
})

test_that("observe_mixture respects weights", {
  set.seed(42)
  obs <- observe_mixture(
    observe_right_censor(tau = 100),
    observe_left_censor(tau = 50),
    weights = c(0.9, 0.1)
  )

  omegas <- character(1000)
  for (i in seq_len(1000)) {
    result <- obs(30)
    omegas[i] <- result$omega
  }

  # With 90% weight on right_censor (which gives "exact" for 30 < 100),
  # most observations should be "exact"
  exact_prop <- mean(omegas == "exact")
  expect_gt(exact_prop, 0.8)
  expect_lt(exact_prop, 0.97)
})

test_that("observe_mixture with uniform weights (default)", {
  set.seed(42)
  obs <- observe_mixture(
    observe_right_censor(tau = 100),
    observe_left_censor(tau = 50)
  )

  omegas <- character(1000)
  for (i in seq_len(1000)) {
    result <- obs(30)
    omegas[i] <- result$omega
  }

  # Roughly 50-50
  exact_prop <- mean(omegas == "exact")
  expect_gt(exact_prop, 0.4)
  expect_lt(exact_prop, 0.6)
})

test_that("observe_mixture errors on empty schemes", {
  expect_error(observe_mixture(), "at least one")
})

test_that("observe_mixture errors on mismatched weights", {
  expect_error(
    observe_mixture(
      observe_right_censor(tau = 10),
      observe_left_censor(tau = 10),
      weights = c(0.5, 0.3, 0.2)
    ),
    "same length"
  )
})

test_that("observe_mixture with three schemes including periodic", {
  set.seed(42)
  obs <- observe_mixture(
    observe_right_censor(tau = 100),
    observe_left_censor(tau = 50),
    observe_periodic(delta = 10, tau = 100),
    weights = c(0.3, 0.3, 0.4)
  )

  omegas <- character(500)
  for (i in seq_len(500)) {
    result <- obs(30)
    omegas[i] <- result$omega
  }

  # Should see exact (from right_censor), left (from left_censor),
  # and interval (from periodic)
  expect_true("exact" %in% omegas)
  expect_true("left" %in% omegas)
  expect_true("interval" %in% omegas)
})


# =============================================================================
# 5. Backwards compatibility: observe=NULL matches old behavior
# =============================================================================

test_that("observe=NULL in generate_masked_series_data matches old behavior", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df_old <- gen(theta = c(0.5, 0.3, 0.2), n = 200, tau = 5, p = 0.3)

  set.seed(42)
  df_new <- gen(
    theta = c(0.5, 0.3, 0.2), n = 200, tau = 5, p = 0.3,
    observe = NULL
  )

  expect_identical(df_old, df_new)
})

test_that("observe=observe_right_censor(tau) matches observe=NULL", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df_null <- gen(theta = c(0.5, 0.3, 0.2), n = 200, tau = 5, p = 0.3)

  set.seed(42)
  df_explicit <- gen(
    theta = c(0.5, 0.3, 0.2), n = 200, tau = 5, p = 0.3,
    observe = observe_right_censor(tau = 5)
  )

  expect_identical(df_null, df_explicit)
})


# =============================================================================
# 6. Integration: rdata with observe functors generates valid data frames
# =============================================================================

test_that("rdata with observe_periodic generates interval data", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 500, p = 0.3,
    observe = observe_periodic(delta = 0.5, tau = 5)
  )

  expect_s3_class(df, "data.frame")
  expect_true("omega" %in% names(df))
  expect_true("t_upper" %in% names(df))
  expect_true(any(df$omega == "interval"))

  # t_upper should be t + delta for interval observations
  interval_rows <- df$omega == "interval"
  if (any(interval_rows)) {
    expect_true(all(df$t_upper[interval_rows] > df$t[interval_rows]))
  }
})

test_that("rdata with observe_left_censor generates left-censored data", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 500, p = 0.3,
    observe = observe_left_censor(tau = 3)
  )

  expect_s3_class(df, "data.frame")
  expect_true(any(df$omega == "left"))
  expect_true(any(df$omega == "right"))
  # No interval observations
  expect_false(any(df$omega == "interval"))
  # All t values should be tau
  expect_true(all(df$t == 3))
})

test_that("rdata with observe_mixture generates mixed-censoring data", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 1000, p = 0.3,
    observe = observe_mixture(
      observe_right_censor(tau = 5),
      observe_left_censor(tau = 3),
      observe_periodic(delta = 0.5, tau = 5),
      weights = c(0.4, 0.3, 0.3)
    )
  )

  expect_s3_class(df, "data.frame")
  omega_types <- unique(df$omega)
  # Should see multiple types
  expect_true(length(omega_types) >= 2)
})

test_that("rdata-generated data with observe is usable by loglik", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 200, p = 0.3,
    observe = observe_periodic(delta = 0.5, tau = 5)
  )

  ll_fn <- loglik(model)
  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))
  expect_true(is.finite(ll))
})

test_that("Weibull rdata with observe_periodic generates valid data", {
  model <- wei_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(1.5, 100, 2.0, 150), n = 200, p = 0.3,
    observe = observe_periodic(delta = 10, tau = 200)
  )

  expect_s3_class(df, "data.frame")
  expect_true(any(df$omega == "interval"))

  ll_fn <- loglik(model)
  ll <- ll_fn(df, par = c(1.5, 100, 2.0, 150))
  expect_true(is.finite(ll))
})

test_that("Homogeneous Weibull rdata with observe_periodic generates valid data", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(1.5, 100, 150), n = 200, p = 0.3,
    observe = observe_periodic(delta = 10, tau = 200)
  )

  expect_s3_class(df, "data.frame")
  expect_true(any(df$omega == "interval"))

  ll_fn <- loglik(model)
  ll <- ll_fn(df, par = c(1.5, 100, 150))
  expect_true(is.finite(ll))
})


# =============================================================================
# 7. Candidate set generation for non-right observations
# =============================================================================

test_that("left-censored observations get candidate sets (C1 satisfied)", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 500, p = 0.3,
    observe = observe_left_censor(tau = 3)
  )

  left_rows <- df$omega == "left"
  if (any(left_rows)) {
    # Each left-censored obs should have at least one component in candidate set
    cand_sums <- df$x1[left_rows] + df$x2[left_rows] + df$x3[left_rows]
    expect_true(all(cand_sums >= 1))
  }
})

test_that("interval-censored observations get candidate sets (C1 satisfied)", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 500, p = 0.3,
    observe = observe_periodic(delta = 0.5, tau = 5)
  )

  interval_rows <- df$omega == "interval"
  if (any(interval_rows)) {
    cand_sums <- df$x1[interval_rows] + df$x2[interval_rows] +
      df$x3[interval_rows]
    expect_true(all(cand_sums >= 1))
  }
})

test_that("right-censored observations have empty candidate sets", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 500, p = 0.3,
    observe = observe_mixture(
      observe_right_censor(tau = 2),
      observe_left_censor(tau = 3),
      weights = c(0.5, 0.5)
    )
  )

  right_rows <- df$omega == "right"
  if (any(right_rows)) {
    cand_sums <- df$x1[right_rows] + df$x2[right_rows] + df$x3[right_rows]
    expect_true(all(cand_sums == 0))
  }
})


# =============================================================================
# 8. p=0 gives singleton candidate sets for non-right observations
# =============================================================================

test_that("p=0 with left-censoring gives singleton candidate sets", {
  model <- exp_series_md_c1_c2_c3()
  gen <- rdata(model)

  set.seed(42)
  df <- gen(
    theta = c(0.5, 0.3, 0.2), n = 500, p = 0,
    observe = observe_left_censor(tau = 3)
  )

  left_rows <- df$omega == "left"
  if (any(left_rows)) {
    cand_sums <- df$x1[left_rows] + df$x2[left_rows] + df$x3[left_rows]
    expect_true(all(cand_sums == 1))
  }
})


# =============================================================================
# 9. Custom column names with interval-censored data
# =============================================================================

test_that("rdata with custom lifetime_upper uses correct column name", {
  model <- exp_series_md_c1_c2_c3(
    lifetime = "time", lifetime_upper = "time_hi",
    omega = "obs_type", candset = "c"
  )
  gen <- rdata(model)

  set.seed(123)
  df <- gen(
    theta = c(0.5, 0.3), n = 200, p = 0.3,
    observe = observe_periodic(delta = 1, tau = 5)
  )

  # Custom column names should be used
  expect_true("time" %in% names(df))
  expect_true("obs_type" %in% names(df))
  expect_true("time_hi" %in% names(df))
  expect_false("t_upper" %in% names(df))

  # The data should be valid for loglik evaluation
  ll_fn <- loglik(model)
  ll_val <- ll_fn(df, par = c(0.5, 0.3))
  expect_true(is.finite(ll_val))
})

