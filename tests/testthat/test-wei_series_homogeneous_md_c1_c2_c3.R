# Tests for wei_series_homogeneous_md_c1_c2_c3 likelihood model

# Helper function to create homogeneous Weibull test data
create_wei_homog_test_data <- function(n = 100, shape = 1.2,
                                        scales = c(100, 150, 200), seed = 42) {
  set.seed(seed)
  m <- length(scales)

  # Generate component lifetimes from Weibull distributions with common shape
  comp_lifetimes <- matrix(NA, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    comp_lifetimes[, j] <- rweibull(n, shape = shape, scale = scales[j])
  }

  # System lifetime is minimum of component lifetimes
  t <- apply(comp_lifetimes, 1, min)

  # Failed component is the one with minimum lifetime
  failed_comp <- apply(comp_lifetimes, 1, which.min)

  # Create candidate sets - include the failed component
  # plus some random other components
  cand_matrix <- matrix(FALSE, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    cand_matrix[i, failed_comp[i]] <- TRUE
    # Randomly include other components with probability 0.3
    for (j in seq_len(m)) {
      if (j != failed_comp[i] && runif(1) < 0.3) {
        cand_matrix[i, j] <- TRUE
      }
    }
  }

  # Build data frame
  df <- data.frame(t = t, delta = TRUE)
  for (j in seq_len(m)) {
    df[[paste0("x", j)]] <- cand_matrix[, j]
  }

  df
}

# Helper to create censored homogeneous Weibull test data
create_wei_homog_censored_data <- function(n = 100, shape = 1.2,
                                            scales = c(100, 150, 200),
                                            tau = 50, seed = 42) {
  df <- create_wei_homog_test_data(n, shape, scales, seed)

  # Apply right censoring
  censored <- df$t > tau
  df$t[censored] <- tau
  df$delta[censored] <- FALSE

  # For censored observations, empty candidate set
  m <- length(scales)
  for (j in seq_len(m)) {
    df[[paste0("x", j)]][censored] <- FALSE
  }

  df
}


# ==============================================================================
# Tests for model constructor
# ==============================================================================

test_that("wei_series_homogeneous_md_c1_c2_c3 constructor creates correct object", {
  model <- wei_series_homogeneous_md_c1_c2_c3()

  expect_s3_class(model, "wei_series_homogeneous_md_c1_c2_c3")
  expect_s3_class(model, "series_md")
  expect_s3_class(model, "likelihood_model")
  expect_equal(model$lifetime, "t")
  expect_equal(model$indicator, "delta")
  expect_equal(model$candset, "x")
})

test_that("wei_series_homogeneous_md_c1_c2_c3 constructor accepts custom column names", {
  model <- wei_series_homogeneous_md_c1_c2_c3(
    lifetime = "time",
    indicator = "censored",
    candset = "cand"
  )

  expect_equal(model$lifetime, "time")
  expect_equal(model$indicator, "censored")
  expect_equal(model$candset, "cand")
})

test_that("wei_series_homogeneous_md_c1_c2_c3 constructor accepts initial parameters", {
  shape <- 1.2
  scales <- c(100, 150, 200)
  model <- wei_series_homogeneous_md_c1_c2_c3(shape = shape, scales = scales)

  expect_equal(model$shape, shape)
  expect_equal(model$scales, scales)
})


# ==============================================================================
# Tests for loglik method
# ==============================================================================

test_that("loglik.wei_series_homogeneous_md_c1_c2_c3 returns a function", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  expect_type(ll_fn, "closure")
})

test_that("loglik returns finite value for valid parameters and data", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_wei_homog_test_data(n = 50)

  # Parameters: (shape, scale1, scale2, scale3)
  par <- c(1.2, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
  expect_type(ll, "double")
})

test_that("loglik returns -Inf for negative shape parameter", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_wei_homog_test_data(n = 50)

  expect_equal(ll_fn(df, par = c(-1.2, 100, 150, 200)), -Inf)
})

test_that("loglik returns -Inf for negative scale parameters", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_wei_homog_test_data(n = 50)

  expect_equal(ll_fn(df, par = c(1.2, -100, 150, 200)), -Inf)
  expect_equal(ll_fn(df, par = c(1.2, 100, -150, 200)), -Inf)
})

test_that("loglik returns -Inf for zero parameters", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_wei_homog_test_data(n = 50)

  expect_equal(ll_fn(df, par = c(0, 100, 150, 200)), -Inf)
  expect_equal(ll_fn(df, par = c(1.2, 0, 150, 200)), -Inf)
})

test_that("loglik errors on empty data frame", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- data.frame(t = numeric(0), x1 = logical(0), x2 = logical(0))

  expect_error(ll_fn(df, par = c(1, 100, 100)), "df is empty")
})

test_that("loglik errors when lifetime column is missing", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- data.frame(x1 = c(TRUE, FALSE), x2 = c(FALSE, TRUE))

  expect_error(ll_fn(df, par = c(1, 100, 100)), "lifetime variable")
})

test_that("loglik errors when parameter length is wrong", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_wei_homog_test_data(n = 50)

  # 3 components need m+1 = 4 parameters, giving only 3
  expect_error(ll_fn(df, par = c(1.2, 100, 150)), "number of components \\+ 1")
})

test_that("loglik handles censored data correctly", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_wei_homog_censored_data(n = 50, tau = 50)

  par <- c(1.2, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("loglik works with backwards compatibility (no delta column)", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # Create data without delta column
  df <- create_wei_homog_test_data(n = 50)
  df$delta <- NULL

  par <- c(1.2, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("loglik value increases when parameters approach true values", {
  # Given: data generated from known parameters
  true_shape <- 1.2
  true_scales <- c(100, 150, 200)
  true_par <- c(1.2, 100, 150, 200)
  distant_par <- c(2.0, 500, 500, 500)

  df <- create_wei_homog_test_data(n = 200, shape = true_shape,
                                    scales = true_scales, seed = 123)

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # When: we evaluate at true params vs distant params
  ll_true <- ll_fn(df, par = true_par)
  ll_distant <- ll_fn(df, par = distant_par)

  # Then: likelihood at true params should be higher
  expect_gt(ll_true, ll_distant)
})


# ==============================================================================
# Tests for score method
# ==============================================================================

test_that("score.wei_series_homogeneous_md_c1_c2_c3 returns a function", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  score_fn <- score(model)

  expect_type(score_fn, "closure")
})

test_that("score returns vector of correct length", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  score_fn <- score(model)
  df <- create_wei_homog_test_data(n = 50)

  # 3 components = 4 parameters (1 shape + 3 scales)
  par <- c(1.2, 100, 150, 200)
  s <- score_fn(df, par = par)

  expect_length(s, 4)
  expect_true(all(is.finite(s)))
})

test_that("score returns NA for negative parameters", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  score_fn <- score(model)
  df <- create_wei_homog_test_data(n = 50)

  s <- score_fn(df, par = c(-1.2, 100, 150, 200))

  expect_true(all(is.na(s)))
})

test_that("score is consistent with numerical gradient", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)

  df <- create_wei_homog_test_data(n = 50, seed = 789)
  par <- c(1.2, 100, 150, 200)

  # Compute analytical score
  analytical <- score_fn(df, par)

  # Compute numerical gradient
  numerical <- numDeriv::grad(func = function(theta) ll_fn(df, theta), x = par)

  # They should be close
  expect_equal(analytical, numerical, tolerance = 1e-4)
})

test_that("score should be approximately zero at MLE", {
  skip_on_cran()  # This test is slow

  # Given: we find the MLE
  true_shape <- 1.2
  true_scales <- c(100, 150, 200)
  set.seed(456)
  df <- create_wei_homog_test_data(n = 300, shape = true_shape, scales = true_scales)

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)

  # Find MLE using optim
  result <- optim(
    par = c(1, 120, 120, 120),
    fn = function(theta) -ll_fn(df, theta),
    method = "L-BFGS-B",
    lower = rep(1e-6, 4)
  )

  mle <- result$par

  # When: we evaluate score at MLE
  s <- score_fn(df, par = mle)

  # Then: score should be approximately zero
  expect_true(all(abs(s) < 0.5),
              info = paste("Score at MLE:", paste(round(s, 4), collapse = ", ")))
})


# ==============================================================================
# Tests for hess_loglik method
# ==============================================================================

test_that("hess_loglik.wei_series_homogeneous_md_c1_c2_c3 returns a function", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)

  expect_type(hess_fn, "closure")
})

test_that("hessian returns matrix of correct dimensions", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)
  df <- create_wei_homog_test_data(n = 50)

  par <- c(1.2, 100, 150, 200)
  H <- hess_fn(df, par = par)

  expect_true(is.matrix(H))
  expect_equal(dim(H), c(4, 4))
  expect_true(all(is.finite(H)))
})

test_that("hessian is approximately symmetric", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)
  df <- create_wei_homog_test_data(n = 50)

  par <- c(1.2, 100, 150, 200)
  H <- hess_fn(df, par = par)

  expect_equal(H, t(H), tolerance = 1e-6)
})

test_that("hessian is consistent with numerical second derivative", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  hess_fn <- hess_loglik(model)

  df <- create_wei_homog_test_data(n = 50, seed = 202)
  par <- c(1.2, 100, 150, 200)

  # Compute analytical Hessian (via jacobian of score)
  analytical <- hess_fn(df, par)

  # Compute numerical Hessian
  numerical <- numDeriv::hessian(func = function(theta) ll_fn(df, theta), x = par)

  # They should be close
  expect_equal(analytical, numerical, tolerance = 1e-3)
})


# ==============================================================================
# Tests for assumptions method
# ==============================================================================

test_that("assumptions returns character vector", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  a <- assumptions(model)

  expect_type(a, "character")
  expect_gt(length(a), 0)
})

test_that("assumptions includes expected assumptions", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  a <- assumptions(model)

  expect_true(any(grepl("Weibull", a, ignore.case = TRUE)))
  expect_true(any(grepl("series", a, ignore.case = TRUE)))
  expect_true(any(grepl("COMMON shape", a, ignore.case = TRUE)))
  expect_true(any(grepl("C1", a)))
  expect_true(any(grepl("C2", a)))
  expect_true(any(grepl("C3", a)))
})


# ==============================================================================
# Tests for wei_series_system_scale helper function
# ==============================================================================

test_that("wei_series_system_scale returns correct system scale", {
  # For k=1 (exponential), system scale = 1 / sum(1/scale_j) = harmonic mean
  k <- 1
  scales <- c(100, 100, 100)

  system_scale <- wei_series_system_scale(k, scales)

  # For k=1, (sum(scale^(-1)))^(-1) = 1/(1/100 + 1/100 + 1/100) = 100/3
  expected <- 100/3
  expect_equal(system_scale, expected, tolerance = 1e-10)
})

test_that("wei_series_system_scale handles different shape values", {
  scales <- c(100, 150, 200)

  # Test with k = 1.5
  k <- 1.5
  system_scale <- wei_series_system_scale(k, scales)

  # Manual calculation
  expected <- (sum(scales^(-k)))^(-1/k)
  expect_equal(system_scale, expected, tolerance = 1e-10)
})

test_that("wei_series_system_scale handles single component", {
  # For a single component system, system scale equals component scale
  k <- 1.2
  scales <- 100

  system_scale <- wei_series_system_scale(k, scales)

  expect_equal(system_scale, 100, tolerance = 1e-10)
})

test_that("wei_series_system_scale decreases with more components", {
  # Adding more components should decrease system lifetime / scale
  k <- 1.2
  scales_2 <- c(100, 150)
  scales_3 <- c(100, 150, 200)
  scales_4 <- c(100, 150, 200, 250)

  scale_2 <- wei_series_system_scale(k, scales_2)
  scale_3 <- wei_series_system_scale(k, scales_3)
  scale_4 <- wei_series_system_scale(k, scales_4)

  expect_gt(scale_2, scale_3)
  expect_gt(scale_3, scale_4)
})


# ==============================================================================
# Tests for edge cases
# ==============================================================================

test_that("loglik handles singleton candidate sets", {
  # Create data where each candidate set contains exactly one component
  df <- data.frame(
    t = c(10, 20, 30),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, FALSE, FALSE),
    x2 = c(FALSE, TRUE, FALSE),
    x3 = c(FALSE, FALSE, TRUE)
  )

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  par <- c(1.2, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("loglik handles all components in candidate set", {
  # Create data where candidate set contains all components
  df <- data.frame(
    t = c(10, 20, 30),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, TRUE, TRUE),
    x2 = c(TRUE, TRUE, TRUE),
    x3 = c(TRUE, TRUE, TRUE)
  )

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  par <- c(1.2, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("model works with custom column names", {
  # Create data with non-standard column names
  df <- data.frame(
    time = c(10, 20, 30),
    censor = c(TRUE, TRUE, TRUE),
    cand1 = c(TRUE, FALSE, FALSE),
    cand2 = c(FALSE, TRUE, TRUE),
    cand3 = c(TRUE, TRUE, TRUE)
  )

  model <- wei_series_homogeneous_md_c1_c2_c3(
    lifetime = "time",
    indicator = "censor",
    candset = "cand"
  )

  ll_fn <- loglik(model)
  par <- c(1.2, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("model handles two-component systems", {
  df <- data.frame(
    t = c(10, 20, 30),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, TRUE, FALSE),
    x2 = c(FALSE, TRUE, TRUE)
  )

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)
  hess_fn <- hess_loglik(model)

  par <- c(1.2, 100, 150)

  ll <- ll_fn(df, par)
  s <- score_fn(df, par)
  H <- hess_fn(df, par)

  expect_true(is.finite(ll))
  expect_length(s, 3)
  expect_equal(dim(H), c(3, 3))
})

test_that("model handles single observation", {
  df <- data.frame(
    t = 15.5,
    delta = TRUE,
    x1 = TRUE,
    x2 = FALSE,
    x3 = TRUE
  )

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  par <- c(1.2, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("model handles shape = 1 (exponential special case)", {
  # When shape is 1, homogeneous Weibull reduces to exponential
  df <- create_wei_homog_test_data(n = 50, shape = 1, scales = c(100, 150, 200))

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  par <- c(1, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("model handles shape < 1 (decreasing hazard)", {
  df <- create_wei_homog_test_data(n = 50, shape = 0.7, scales = c(100, 150, 200))

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  par <- c(0.7, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})

test_that("model handles shape > 1 (increasing hazard)", {
  df <- create_wei_homog_test_data(n = 50, shape = 2.5, scales = c(100, 150, 200))

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  par <- c(2.5, 100, 150, 200)
  ll <- ll_fn(df, par = par)

  expect_true(is.finite(ll))
})


# ==============================================================================
# Tests comparing homogeneous to full Weibull model
# ==============================================================================

test_that("homogeneous model gives same likelihood as full model when shapes are equal", {
  # When all shapes are the same, homogeneous and full models should agree
  shape <- 1.2
  scales <- c(100, 150, 200)

  df <- create_wei_homog_test_data(n = 50, shape = shape, scales = scales, seed = 999)

  # Homogeneous model
  model_homog <- wei_series_homogeneous_md_c1_c2_c3()
  ll_homog <- loglik(model_homog)
  par_homog <- c(shape, scales)

  # Full model (with same shape for all components)
  model_full <- wei_series_md_c1_c2_c3()
  ll_full <- loglik(model_full)
  par_full <- c(shape, scales[1], shape, scales[2], shape, scales[3])

  ll_h <- ll_homog(df, par_homog)
  ll_f <- ll_full(df, par_full)

  expect_equal(ll_h, ll_f, tolerance = 1e-10)
})


# ==============================================================================
# Tests for rdata method
# ==============================================================================

test_that("rdata.wei_series_homogeneous_md_c1_c2_c3 returns a function", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  expect_type(rdata_fn, "closure")
})

test_that("rdata generates data frame with correct structure", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  set.seed(42)
  # theta = (k, scale1, scale2, scale3)
  df <- rdata_fn(theta = c(1.5, 100, 150, 200), n = 100)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 100)
  expect_true("t" %in% names(df))
  expect_true("delta" %in% names(df))
  expect_true("x1" %in% names(df))
  expect_true("x2" %in% names(df))
  expect_true("x3" %in% names(df))
})

test_that("rdata respects custom column names from model", {
  model <- wei_series_homogeneous_md_c1_c2_c3(
    lifetime = "time",
    indicator = "censored",
    candset = "cand"
  )
  rdata_fn <- rdata(model)

  set.seed(42)
  # theta = (k, scale1, scale2)
  df <- rdata_fn(theta = c(1.5, 100, 150), n = 50)

  expect_true("time" %in% names(df))
  expect_true("censored" %in% names(df))
  expect_true("cand1" %in% names(df))
  expect_true("cand2" %in% names(df))
})

test_that("rdata applies right-censoring correctly", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  set.seed(42)
  df <- rdata_fn(theta = c(1.5, 100, 150, 200), n = 200, tau = 30)

  # Check that censored observations exist and have correct properties
  censored <- !df$delta
  if (any(censored)) {
    # Censored times should equal tau
    expect_true(all(df$t[censored] == 30))
    # Censored observations should have empty candidate sets
    expect_true(all(!df$x1[censored]))
    expect_true(all(!df$x2[censored]))
    expect_true(all(!df$x3[censored]))
  }
})

test_that("rdata generates candidate sets satisfying C1", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  # With p = 0, only the failed component should be in candidate set
  set.seed(42)
  df <- rdata_fn(theta = c(1.5, 100, 150, 200), n = 100, p = 0)

  # Each exact observation should have exactly one component in candidate set
  exact <- df$delta
  cand_sums <- df$x1[exact] + df$x2[exact] + df$x3[exact]
  expect_true(all(cand_sums == 1))
})

test_that("rdata errors on invalid parameters", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  # Negative shape
  expect_error(rdata_fn(theta = c(-1, 100, 150), n = 100))
  # Zero scale
  expect_error(rdata_fn(theta = c(1.5, 0, 150), n = 100))
})

test_that("rdata generated data can be used with loglik", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rdata_fn <- rdata(model)
  ll_fn <- loglik(model)

  set.seed(42)
  theta <- c(1.5, 100, 150, 200)
  df <- rdata_fn(theta = theta, n = 100, tau = 300, p = 0.3)
  ll <- ll_fn(df, par = theta)

  expect_true(is.finite(ll))
})


# ==============================================================================
# Tests for observed_info method
# ==============================================================================

test_that("observed_info returns a function", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  obs_info_fn <- observed_info(model)

  expect_type(obs_info_fn, "closure")
})

test_that("observed_info returns correct dimensions", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  obs_info_fn <- observed_info(model)

  df <- create_wei_homog_test_data(n = 50)
  par <- c(1.5, 100, 150, 200)

  I_obs <- obs_info_fn(df, par)

  expect_true(is.matrix(I_obs))
  expect_equal(dim(I_obs), c(4, 4))
})


# ==============================================================================
# Tests for fim method
# ==============================================================================

test_that("fim returns a function", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  fim_fn <- fim(model)

  expect_type(fim_fn, "closure")
})

test_that("fim returns matrix of correct dimensions", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  fim_fn <- fim(model)

  set.seed(42)
  theta <- c(1.5, 100, 150, 200)
  I <- fim_fn(theta = theta, n_obs = 50, n_samples = 30)

  expect_true(is.matrix(I))
  expect_equal(dim(I), c(4, 4))
  expect_true(all(is.finite(I)))
})


# ==============================================================================
# Test: Single-component system (m=1)
# ==============================================================================

test_that("homogeneous Weibull works with single-component system (m=1)", {
  df <- data.frame(
    t = c(10, 20, 30),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, TRUE, TRUE)
  )

  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)
  hess_fn <- hess_loglik(model)

  par <- c(1.2, 100)

  ll <- ll_fn(df, par = par)
  s <- score_fn(df, par = par)
  H <- hess_fn(df, par = par)

  expect_true(is.finite(ll))
  expect_length(s, 2)
  expect_equal(dim(H), c(2, 2))
})


# ==============================================================================
# Test: High censoring fraction
# ==============================================================================

test_that("homogeneous Weibull loglik is finite with high censoring", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  set.seed(42)
  df <- rdata_fn(theta = c(1.5, 100, 150, 200), n = 200, tau = 5, p = 0.3)
  cens_frac <- mean(!df$delta)
  expect_gt(cens_frac, 0.5)

  ll_fn <- loglik(model)
  score_fn <- score(model)

  par <- c(1.5, 100, 150, 200)
  ll <- ll_fn(df, par = par)
  s <- score_fn(df, par = par)

  expect_true(is.finite(ll))
  expect_true(all(is.finite(s)))
})


# ==============================================================================
# Test: wei_series_system_scale edge cases
# ==============================================================================

test_that("wei_series_system_scale with k=1 gives harmonic-mean-like result", {
  scales <- c(100, 200, 300)
  result <- wei_series_system_scale(k = 1, scales = scales)
  expected <- 1 / sum(1 / scales)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("wei_series_system_scale with single scale returns that scale", {
  k <- 2.5
  scale <- 150
  result <- wei_series_system_scale(k, scale)
  expect_equal(result, scale, tolerance = 1e-10)
})

test_that("wei_series_system_scale handles very large k", {
  # As k -> Inf, system scale -> min(scales)
  scales <- c(100, 150, 200)
  result <- wei_series_system_scale(k = 100, scales = scales)
  expect_equal(result, min(scales), tolerance = 1)
})

test_that("wei_series_system_scale handles very small k", {
  # As k -> 0+, behavior should still be finite
  scales <- c(100, 150, 200)
  result <- wei_series_system_scale(k = 0.01, scales = scales)
  expect_true(is.finite(result))
  expect_true(result > 0)
})


# ==============================================================================
# Test: Error message validation
# ==============================================================================

test_that("homogeneous Weibull error messages are informative", {
  model <- wei_series_homogeneous_md_c1_c2_c3()
  ll_fn <- loglik(model)

  df_empty <- data.frame(t = numeric(0), x1 = logical(0))
  expect_error(ll_fn(df_empty, par = c(1, 100)), "df is empty")

  df_no_t <- data.frame(x1 = c(TRUE, FALSE))
  expect_error(ll_fn(df_no_t, par = c(1, 100)), "lifetime variable")

  df <- data.frame(t = c(1, 2), delta = c(TRUE, TRUE),
                   x1 = c(TRUE, FALSE), x2 = c(FALSE, TRUE))
  expect_error(ll_fn(df, par = c(1.2, 100)), "number of components \\+ 1")
})
