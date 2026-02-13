# Tests for series_md generics:
#   ncomponents, component_hazard, conditional_cause_probability, cause_probability

# =============================================================================
# ncomponents
# =============================================================================

test_that("ncomponents returns correct m for exponential model with rates", {
  model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  expect_equal(ncomponents(model), 3L)
})

test_that("ncomponents returns NULL for exponential model without rates", {
  model <- exp_series_md_c1_c2_c3()
  expect_null(ncomponents(model))
})

test_that("ncomponents returns correct m for Weibull model with shapes", {
  model <- wei_series_md_c1_c2_c3(shapes = c(1.5, 2.0), scales = c(100, 200))
  expect_equal(ncomponents(model), 2L)
})

test_that("ncomponents returns NULL for Weibull model without shapes", {
  model <- wei_series_md_c1_c2_c3()
  expect_null(ncomponents(model))
})

test_that("ncomponents returns correct m for homogeneous Weibull with scales", {
  model <- wei_series_homogeneous_md_c1_c2_c3(
    shape = 1.5, scales = c(100, 200, 300)
  )
  expect_equal(ncomponents(model), 3L)
})

test_that("ncomponents returns NULL for homogeneous Weibull without scales", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  expect_null(ncomponents(model))
})


# =============================================================================
# component_hazard — exponential model
# =============================================================================

test_that("component_hazard returns closure for exponential model", {
  model <- exp_series_md_c1_c2_c3()
  h <- component_hazard(model, 1)
  expect_type(h, "closure")
})

test_that("exponential component_hazard returns constant rate (memoryless)", {
  model <- exp_series_md_c1_c2_c3()
  rates <- c(0.5, 0.3, 0.2)
  t_vals <- c(1, 5, 10, 100)

  for (j in seq_along(rates)) {
    h_j <- component_hazard(model, j)
    result <- h_j(t_vals, par = rates)
    expect_equal(result, rep(rates[j], length(t_vals)))
  }
})


# =============================================================================
# component_hazard — Weibull model
# =============================================================================

test_that("component_hazard returns closure for Weibull model", {
  model <- wei_series_md_c1_c2_c3()
  h <- component_hazard(model, 1)
  expect_type(h, "closure")
})

test_that("Weibull component_hazard matches hand-computed values", {
  model <- wei_series_md_c1_c2_c3()
  # par = (k1, beta1, k2, beta2) = (2, 100, 1.5, 200)
  par <- c(2, 100, 1.5, 200)
  t_val <- 50

  # h_1(50) = (2/100) * (50/100)^(2-1) = 0.02 * 0.5 = 0.01
  h1 <- component_hazard(model, 1)
  expect_equal(h1(t_val, par), 0.01)

  # h_2(50) = (1.5/200) * (50/200)^(1.5-1) = 0.0075 * 0.25^0.5 = 0.0075 * 0.5
  h2 <- component_hazard(model, 2)
  expect_equal(h2(t_val, par), 0.0075 * 0.5, tolerance = 1e-10)
})

test_that("Weibull(shape=1) hazard equals exponential rate", {
  model <- wei_series_md_c1_c2_c3()
  rate <- 0.5
  scale <- 1 / rate
  par <- c(1, scale)
  t_vals <- c(1, 5, 10)

  h <- component_hazard(model, 1)
  result <- h(t_vals, par)
  expect_equal(result, rep(rate, 3), tolerance = 1e-10)
})


# =============================================================================
# component_hazard — homogeneous Weibull model
# =============================================================================

test_that("component_hazard returns closure for homogeneous Weibull model", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  h <- component_hazard(model, 1)
  expect_type(h, "closure")
})

test_that("homogeneous Weibull component_hazard matches hand-computed", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  # par = (k, beta1, beta2) = (2, 100, 200)
  par <- c(2, 100, 200)
  t_val <- 50

  # h_1(50) = (2/100) * (50/100)^(2-1) = 0.02 * 0.5 = 0.01
  h1 <- component_hazard(model, 1)
  expect_equal(h1(t_val, par), 0.01)

  # h_2(50) = (2/200) * (50/200)^(2-1) = 0.01 * 0.25 = 0.0025
  h2 <- component_hazard(model, 2)
  expect_equal(h2(t_val, par), 0.0025)
})

test_that("homogeneous Weibull(shape=1) hazard equals exponential rate", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rates <- c(0.5, 0.3)
  scales <- 1 / rates
  par <- c(1, scales)

  for (j in seq_along(rates)) {
    h_j <- component_hazard(model, j)
    result <- h_j(c(1, 5, 10), par)
    expect_equal(result, rep(rates[j], 3), tolerance = 1e-10)
  }
})


# =============================================================================
# conditional_cause_probability — general properties
# =============================================================================

test_that("conditional_cause_probability returns closure for each model", {
  exp_model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  wei_model <- wei_series_md_c1_c2_c3(shapes = c(1.5, 2), scales = c(100, 200))
  hom_model <- wei_series_homogeneous_md_c1_c2_c3(
    shape = 1.5, scales = c(100, 200, 300)
  )

  expect_type(conditional_cause_probability(exp_model), "closure")
  expect_type(conditional_cause_probability(wei_model), "closure")
  expect_type(conditional_cause_probability(hom_model), "closure")
})

test_that("conditional probabilities are in [0,1] and rows sum to 1", {
  models_and_pars <- list(
    list(
      model = exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2)),
      par = c(0.5, 0.3, 0.2),
      m = 3
    ),
    list(
      model = wei_series_md_c1_c2_c3(shapes = c(1.5, 2), scales = c(100, 200)),
      par = c(1.5, 100, 2, 200),
      m = 2
    ),
    list(
      model = wei_series_homogeneous_md_c1_c2_c3(
        shape = 1.5, scales = c(100, 200)
      ),
      par = c(1.5, 100, 200),
      m = 2
    )
  )

  t_vals <- c(1, 10, 50, 100)

  for (spec in models_and_pars) {
    ccp_fn <- conditional_cause_probability(spec$model)
    probs <- ccp_fn(t_vals, spec$par)

    expect_equal(nrow(probs), length(t_vals))
    expect_equal(ncol(probs), spec$m)
    expect_true(all(probs >= 0))
    expect_true(all(probs <= 1))
    expect_equal(rowSums(probs), rep(1, length(t_vals)), tolerance = 1e-10)
  }
})


# =============================================================================
# conditional_cause_probability — exponential-specific
# =============================================================================

test_that("exponential conditional cause probability is time-independent", {
  model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  ccp_fn <- conditional_cause_probability(model)
  par <- c(0.5, 0.3, 0.2)

  probs <- ccp_fn(c(0.01, 1, 10, 100, 1000), par)

  # All rows should be identical
  for (i in 2:nrow(probs)) {
    expect_equal(probs[i, ], probs[1, ], tolerance = 1e-15)
  }
})

test_that("exponential conditional probability equals lambda_j / sum(lambda)", {
  model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  ccp_fn <- conditional_cause_probability(model)
  par <- c(0.5, 0.3, 0.2)

  probs <- ccp_fn(5.0, par)
  expected <- par / sum(par)

  expect_equal(as.numeric(probs), expected, tolerance = 1e-15)
})


# =============================================================================
# conditional_cause_probability — Weibull varies with t
# =============================================================================

test_that("Weibull conditional cause probability varies with t", {
  model <- wei_series_md_c1_c2_c3(shapes = c(0.5, 2.0), scales = c(100, 100))
  ccp_fn <- conditional_cause_probability(model)
  par <- c(0.5, 100, 2.0, 100)

  probs_early <- ccp_fn(1, par)
  probs_late <- ccp_fn(200, par)

  # With shape 0.5 (decreasing hazard) vs shape 2.0 (increasing hazard),
  # early failures favor component 1, late failures favor component 2
  expect_gt(probs_early[1, 1], probs_late[1, 1])
  expect_lt(probs_early[1, 2], probs_late[1, 2])
})


# =============================================================================
# conditional_cause_probability — Weibull hand-computed
# =============================================================================

test_that("Weibull conditional probability matches hand-computed values", {
  model <- wei_series_md_c1_c2_c3(shapes = c(2, 2), scales = c(100, 200))
  ccp_fn <- conditional_cause_probability(model)
  par <- c(2, 100, 2, 200)
  t_val <- 50

  # h_1(50) = (2/100) * (50/100)^1 = 0.01
  # h_2(50) = (2/200) * (50/200)^1 = 0.0025
  # P(K=1|T=50) = 0.01 / 0.0125 = 0.8
  # P(K=2|T=50) = 0.0025 / 0.0125 = 0.2
  probs <- ccp_fn(t_val, par)
  expect_equal(probs[1, 1], 0.8, tolerance = 1e-10)
  expect_equal(probs[1, 2], 0.2, tolerance = 1e-10)
})


# =============================================================================
# cause_probability — general properties
# =============================================================================

test_that("cause_probability returns closure for each model", {
  exp_model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  wei_model <- wei_series_md_c1_c2_c3(shapes = c(1.5, 2), scales = c(100, 200))
  hom_model <- wei_series_homogeneous_md_c1_c2_c3(
    shape = 1.5, scales = c(100, 200)
  )

  expect_type(cause_probability(exp_model), "closure")
  expect_type(cause_probability(wei_model), "closure")
  expect_type(cause_probability(hom_model), "closure")
})

test_that("cause probabilities are in [0,1] and sum to 1", {
  models_and_pars <- list(
    list(
      model = exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2)),
      par = c(0.5, 0.3, 0.2),
      m = 3
    ),
    list(
      model = wei_series_md_c1_c2_c3(shapes = c(1.5, 2), scales = c(100, 200)),
      par = c(1.5, 100, 2, 200),
      m = 2
    ),
    list(
      model = wei_series_homogeneous_md_c1_c2_c3(
        shape = 1.5, scales = c(100, 200)
      ),
      par = c(1.5, 100, 200),
      m = 2
    )
  )

  for (spec in models_and_pars) {
    set.seed(42)
    cp_fn <- cause_probability(spec$model)
    probs <- cp_fn(spec$par, n_mc = 5000)

    expect_length(probs, spec$m)
    expect_true(all(probs >= 0))
    expect_true(all(probs <= 1))
    expect_equal(sum(probs), 1, tolerance = 0.02)
  }
})


# =============================================================================
# cause_probability — exponential-specific
# =============================================================================

test_that("exponential cause_probability equals lambda_j / sum(lambda)", {
  model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  cp_fn <- cause_probability(model)
  par <- c(0.5, 0.3, 0.2)

  probs <- cp_fn(par)
  expected <- par / sum(par)

  expect_equal(probs, expected, tolerance = 1e-15)
})

test_that("exponential cause_probability matches conditional (time-independent)", {
  model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  par <- c(0.5, 0.3, 0.2)

  cp_fn <- cause_probability(model)
  ccp_fn <- conditional_cause_probability(model)

  cp <- cp_fn(par)
  ccp <- ccp_fn(42.0, par)  # arbitrary t

  expect_equal(cp, as.numeric(ccp), tolerance = 1e-15)
})


# =============================================================================
# cause_probability — Weibull Monte Carlo
# =============================================================================

test_that("Weibull cause_probability MC estimate is reasonable", {
  skip_on_cran()

  model <- wei_series_md_c1_c2_c3(shapes = c(2, 2), scales = c(100, 200))
  cp_fn <- cause_probability(model)
  par <- c(2, 100, 2, 200)

  # When shapes are equal, P(K=j) = beta_j^{-k} / sum(beta_l^{-k})
  # = 100^{-2} / (100^{-2} + 200^{-2}) = 1e-4 / 1.25e-4 = 0.8
  expected <- c(100^(-2), 200^(-2))
  expected <- expected / sum(expected)

  set.seed(42)
  probs <- cp_fn(par, n_mc = 50000)

  expect_equal(probs, expected, tolerance = 0.02)
})

test_that("homogeneous Weibull cause_probability MC estimate is reasonable", {
  skip_on_cran()

  model <- wei_series_homogeneous_md_c1_c2_c3(
    shape = 2, scales = c(100, 200, 300)
  )
  cp_fn <- cause_probability(model)
  par <- c(2, 100, 200, 300)

  # P(K=j) = beta_j^{-k} / sum(beta_l^{-k})
  inv_scales_k <- c(100, 200, 300)^(-2)
  expected <- inv_scales_k / sum(inv_scales_k)

  set.seed(42)
  probs <- cp_fn(par, n_mc = 50000)

  expect_equal(probs, expected, tolerance = 0.02)
})


# =============================================================================
# Consistency: averaged conditional ≈ cause_probability
# =============================================================================

test_that("averaged conditional ≈ cause_probability for exponential", {
  model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  par <- c(0.5, 0.3, 0.2)

  ccp_fn <- conditional_cause_probability(model)
  cp_fn <- cause_probability(model)

  # Generate many system lifetimes
  set.seed(42)
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = par, n = 5000)
  t_vals <- df$t[df$omega == "exact"]

  avg_cond <- colMeans(ccp_fn(t_vals, par))
  marginal <- cp_fn(par)

  expect_equal(avg_cond, marginal, tolerance = 1e-10)
})

test_that("averaged conditional ≈ cause_probability for Weibull", {
  skip_on_cran()

  model <- wei_series_md_c1_c2_c3(shapes = c(1.5, 2), scales = c(100, 200))
  par <- c(1.5, 100, 2, 200)

  ccp_fn <- conditional_cause_probability(model)
  cp_fn <- cause_probability(model)

  set.seed(42)
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = par, n = 50000)
  t_vals <- df$t[df$omega == "exact"]

  avg_cond <- colMeans(ccp_fn(t_vals, par))

  set.seed(123)
  marginal <- cp_fn(par, n_mc = 50000)

  expect_equal(avg_cond, marginal, tolerance = 0.02)
})


# =============================================================================
# Cross-model: Weibull(shape=1) matches exponential for all generics
# =============================================================================

test_that("Weibull(shape=1) component_hazard matches exponential", {
  rates <- c(0.5, 0.3)
  scales <- 1 / rates
  t_vals <- c(1, 5, 10)

  exp_model <- exp_series_md_c1_c2_c3()
  wei_model <- wei_series_md_c1_c2_c3()
  hom_model <- wei_series_homogeneous_md_c1_c2_c3()

  exp_par <- rates
  wei_par <- c(1, scales[1], 1, scales[2])
  hom_par <- c(1, scales)

  for (j in seq_along(rates)) {
    h_exp <- component_hazard(exp_model, j)(t_vals, exp_par)
    h_wei <- component_hazard(wei_model, j)(t_vals, wei_par)
    h_hom <- component_hazard(hom_model, j)(t_vals, hom_par)

    expect_equal(h_exp, h_wei, tolerance = 1e-10)
    expect_equal(h_exp, h_hom, tolerance = 1e-10)
  }
})

test_that("Weibull(shape=1) conditional cause prob matches exponential", {
  rates <- c(0.5, 0.3)
  scales <- 1 / rates
  t_vals <- c(1, 5, 10)

  exp_model <- exp_series_md_c1_c2_c3(rates = rates)
  wei_model <- wei_series_md_c1_c2_c3(shapes = c(1, 1), scales = scales)
  hom_model <- wei_series_homogeneous_md_c1_c2_c3(shape = 1, scales = scales)

  exp_par <- rates
  wei_par <- c(1, scales[1], 1, scales[2])
  hom_par <- c(1, scales)

  p_exp <- conditional_cause_probability(exp_model)(t_vals, exp_par)
  p_wei <- conditional_cause_probability(wei_model)(t_vals, wei_par)
  p_hom <- conditional_cause_probability(hom_model)(t_vals, hom_par)

  expect_equal(p_exp, p_wei, tolerance = 1e-10)
  expect_equal(p_exp, p_hom, tolerance = 1e-10)
})


# =============================================================================
# Edge cases
# =============================================================================

test_that("component_hazard works with single-component system", {
  exp_model <- exp_series_md_c1_c2_c3()
  h <- component_hazard(exp_model, 1)
  expect_equal(h(5, par = 0.3), 0.3)

  wei_model <- wei_series_md_c1_c2_c3()
  h <- component_hazard(wei_model, 1)
  # h_1(5) = (2/10) * (5/10)^1 = 0.1
  expect_equal(h(5, par = c(2, 10)), 0.1)
})

test_that("conditional_cause_probability handles scalar t", {
  model <- exp_series_md_c1_c2_c3(rates = c(0.5, 0.3, 0.2))
  ccp_fn <- conditional_cause_probability(model)
  par <- c(0.5, 0.3, 0.2)

  probs <- ccp_fn(5.0, par)
  expect_equal(nrow(probs), 1)
  expect_equal(ncol(probs), 3)
  expect_equal(sum(probs), 1, tolerance = 1e-15)
})

test_that("conditional_cause_probability errors for NULL ncomponents", {
  model <- wei_series_md_c1_c2_c3()
  # Default series_md method will be dispatched since no specific override
  # exists for wei_series_md_c1_c2_c3. ncomponents returns NULL, so
  # the closure should error when called.
  ccp_fn <- conditional_cause_probability(model)
  expect_error(ccp_fn(5.0, par = c(1.5, 100, 2, 200)), "ncomponents")
})
