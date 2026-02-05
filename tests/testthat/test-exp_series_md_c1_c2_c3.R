# Tests for exp_series_md_c1_c2_c3 likelihood model

# Helper function to create test data
create_exp_test_data <- function(n = 100, rates = c(0.5, 0.3, 0.2), seed = 42) {
  set.seed(seed)
  m <- length(rates)

  # Generate component lifetimes
  comp_lifetimes <- matrix(NA, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    comp_lifetimes[, j] <- rexp(n, rates[j])
  }

  # System lifetime is minimum of component lifetimes
  t <- apply(comp_lifetimes, 1, min)

  # Failed component is the one with minimum lifetime
  failed_comp <- apply(comp_lifetimes, 1, which.min)

  # Create candidate sets - for simplicity, include the failed component
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

# Helper to create censored test data
create_exp_censored_data <- function(n = 100, rates = c(0.5, 0.3, 0.2),
                                      tau = 2, seed = 42) {
  df <- create_exp_test_data(n, rates, seed)

  # Apply right censoring
  censored <- df$t > tau
  df$t[censored] <- tau
  df$delta[censored] <- FALSE

  # For censored observations, empty candidate set
  m <- length(rates)
  for (j in seq_len(m)) {
    df[[paste0("x", j)]][censored] <- FALSE
  }

  df
}


# ==============================================================================
# Tests for model constructor
# ==============================================================================

test_that("exp_series_md_c1_c2_c3 constructor creates correct object", {
  model <- exp_series_md_c1_c2_c3()

  expect_s3_class(model, "exp_series_md_c1_c2_c3")
  expect_s3_class(model, "series_md")
  expect_s3_class(model, "likelihood_model")
  expect_equal(model$lifetime, "t")
  expect_equal(model$indicator, "delta")
  expect_equal(model$candset, "x")
})

test_that("exp_series_md_c1_c2_c3 constructor accepts custom column names", {
  model <- exp_series_md_c1_c2_c3(
    lifetime = "time",
    indicator = "censored",
    candset = "cand"
  )

  expect_equal(model$lifetime, "time")
  expect_equal(model$indicator, "censored")
  expect_equal(model$candset, "cand")
})

test_that("exp_series_md_c1_c2_c3 constructor accepts initial rates", {
  rates <- c(0.5, 0.3, 0.2)
  model <- exp_series_md_c1_c2_c3(rates = rates)

  expect_equal(model$rates, rates)
})


# ==============================================================================
# Tests for loglik method
# ==============================================================================

test_that("loglik.exp_series_md_c1_c2_c3 returns a function", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  expect_type(ll_fn, "closure")
})

test_that("loglik returns finite value for valid parameters and data", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_exp_test_data(n = 50)

  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
  expect_type(ll, "double")
})

test_that("loglik returns -Inf for negative parameters", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_exp_test_data(n = 50)

  expect_equal(ll_fn(df, par = c(-0.5, 0.3, 0.2)), -Inf)
  expect_equal(ll_fn(df, par = c(0.5, -0.3, 0.2)), -Inf)
  expect_equal(ll_fn(df, par = c(0.5, 0.3, -0.2)), -Inf)
})

test_that("loglik returns -Inf for zero parameters", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_exp_test_data(n = 50)

  expect_equal(ll_fn(df, par = c(0, 0.3, 0.2)), -Inf)
})

test_that("loglik errors on empty data frame", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- data.frame(t = numeric(0), x1 = logical(0), x2 = logical(0))

  expect_error(ll_fn(df, par = c(0.5, 0.3)), "df is empty")
})

test_that("loglik errors when lifetime column is missing", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- data.frame(x1 = c(TRUE, FALSE), x2 = c(FALSE, TRUE))

  expect_error(ll_fn(df, par = c(0.5, 0.3)))
})

test_that("loglik errors when no candidate set columns found", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- data.frame(t = c(1, 2, 3))

  # When no candidate set columns are found, md_decode_matrix returns NULL
  # which causes an error in ncol(C)
  expect_error(ll_fn(df, par = c(0.5)))
})

test_that("loglik handles censored data correctly", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  df <- create_exp_censored_data(n = 50, tau = 1)

  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
})

test_that("loglik works with backwards compatibility (no delta column)", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # Create data without delta column
  df <- create_exp_test_data(n = 50)
  df$delta <- NULL

  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
})

test_that("loglik value increases when parameters approach true values", {
  # Given: data generated from known parameters
  true_rates <- c(0.5, 0.3, 0.2)
  df <- create_exp_test_data(n = 200, rates = true_rates, seed = 123)

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # When: we evaluate at true params vs distant params
  ll_true <- ll_fn(df, par = true_rates)
  ll_distant <- ll_fn(df, par = c(2.0, 2.0, 2.0))


  # Then: likelihood at true params should be higher
  expect_gt(ll_true, ll_distant)
})


# ==============================================================================
# Tests for score method
# ==============================================================================

test_that("score.exp_series_md_c1_c2_c3 returns a function", {
  model <- exp_series_md_c1_c2_c3()
  score_fn <- score(model)

  expect_type(score_fn, "closure")
})

test_that("score returns vector of correct length", {
  model <- exp_series_md_c1_c2_c3()
  score_fn <- score(model)
  df <- create_exp_test_data(n = 50)

  s <- score_fn(df, par = c(0.5, 0.3, 0.2))

  expect_length(s, 3)
  expect_true(all(is.finite(s)))
})

test_that("score returns NA for negative parameters", {
  model <- exp_series_md_c1_c2_c3()
  score_fn <- score(model)
  df <- create_exp_test_data(n = 50)

  s <- score_fn(df, par = c(-0.5, 0.3, 0.2))

  expect_true(all(is.na(s)))
})

test_that("score should be approximately zero at MLE", {
  # Given: we find the MLE
  true_rates <- c(0.5, 0.3, 0.2)
  set.seed(456)
  df <- create_exp_test_data(n = 500, rates = true_rates)

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)

  # Find MLE using optim
  result <- optim(
    par = c(1, 1, 1),
    fn = function(theta) -ll_fn(df, theta),
    method = "L-BFGS-B",
    lower = rep(1e-6, 3)
  )

  mle <- result$par

  # When: we evaluate score at MLE
  s <- score_fn(df, par = mle)

  # Then: score should be approximately zero
  expect_true(all(abs(s) < 0.1),
              info = paste("Score at MLE:", paste(round(s, 4), collapse = ", ")))
})

test_that("score is consistent with numerical gradient", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)

  df <- create_exp_test_data(n = 50, seed = 789)
  par <- c(0.5, 0.3, 0.2)

  # Compute analytical score
  analytical <- score_fn(df, par)

  # Compute numerical gradient
  numerical <- numDeriv::grad(func = function(theta) ll_fn(df, theta), x = par)

  # They should be close
  expect_equal(analytical, numerical, tolerance = 1e-5)
})


# ==============================================================================
# Tests for hess_loglik method
# ==============================================================================

test_that("hess_loglik.exp_series_md_c1_c2_c3 returns a function", {
  model <- exp_series_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)

  expect_type(hess_fn, "closure")
})

test_that("hessian returns matrix of correct dimensions", {
  model <- exp_series_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)
  df <- create_exp_test_data(n = 50)

  H <- hess_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.matrix(H))
  expect_equal(dim(H), c(3, 3))
  expect_true(all(is.finite(H)))
})

test_that("hessian is symmetric", {
  model <- exp_series_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)
  df <- create_exp_test_data(n = 50)

  H <- hess_fn(df, par = c(0.5, 0.3, 0.2))

  expect_equal(H, t(H), tolerance = 1e-10)
})

test_that("hessian should be negative semi-definite at MLE", {
  # Given: we find the MLE
  true_rates <- c(0.5, 0.3, 0.2)
  set.seed(101)
  df <- create_exp_test_data(n = 200, rates = true_rates)

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  hess_fn <- hess_loglik(model)

  # Find MLE
  result <- optim(
    par = c(1, 1, 1),
    fn = function(theta) -ll_fn(df, theta),
    method = "L-BFGS-B",
    lower = rep(1e-6, 3)
  )
  mle <- result$par

  # When: we compute Hessian at MLE
  H <- hess_fn(df, par = mle)

  # Then: all eigenvalues should be non-positive (negative semi-definite)
  eigenvalues <- eigen(H, symmetric = TRUE)$values
  expect_true(all(eigenvalues <= 1e-6),
              info = paste("Eigenvalues:", paste(round(eigenvalues, 4), collapse = ", ")))
})

test_that("hessian returns NA matrix for invalid parameters", {
  model <- exp_series_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)
  df <- create_exp_test_data(n = 50)

  H <- hess_fn(df, par = c(-0.5, 0.3, 0.2))

  expect_true(all(is.na(H)))
})

test_that("hessian is consistent with numerical second derivative", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  hess_fn <- hess_loglik(model)

  df <- create_exp_test_data(n = 50, seed = 202)
  par <- c(0.5, 0.3, 0.2)

  # Compute analytical Hessian
  analytical <- hess_fn(df, par)

  # Compute numerical Hessian
  numerical <- numDeriv::hessian(func = function(theta) ll_fn(df, theta), x = par)

  # They should be close
  expect_equal(analytical, numerical, tolerance = 1e-4)
})


# ==============================================================================
# Tests for assumptions method
# ==============================================================================

test_that("assumptions returns character vector", {
  model <- exp_series_md_c1_c2_c3()
  a <- assumptions(model)

  expect_type(a, "character")
  expect_gt(length(a), 0)
})

test_that("assumptions includes expected assumptions", {
  model <- exp_series_md_c1_c2_c3()
  a <- assumptions(model)

  expect_true(any(grepl("exponential", a, ignore.case = TRUE)))
  expect_true(any(grepl("series", a, ignore.case = TRUE)))
  expect_true(any(grepl("C1", a)))
  expect_true(any(grepl("C2", a)))
  expect_true(any(grepl("C3", a)))
})


# ==============================================================================
# Tests for MLE convergence and consistency
# ==============================================================================

test_that("MLE estimates converge to true parameters with large samples", {
  skip_on_cran()  # This test is slow

  # Given: large sample from known parameters
  true_rates <- c(0.5, 0.3, 0.2)
  set.seed(303)
  df <- create_exp_test_data(n = 1000, rates = true_rates)

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # When: we find MLE
  result <- optim(
    par = c(1, 1, 1),
    fn = function(theta) -ll_fn(df, theta),
    method = "L-BFGS-B",
    lower = rep(1e-6, 3)
  )
  mle <- result$par

  # Then: MLE should be close to true parameters
  # (allow some tolerance due to masking)
  expect_equal(mle, true_rates, tolerance = 0.15)
})


# ==============================================================================
# Tests for edge cases
# ==============================================================================

test_that("loglik handles singleton candidate sets (exact data)", {
  # Create data where each candidate set contains exactly one component
  df <- data.frame(
    t = c(1, 2, 3),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, FALSE, FALSE),
    x2 = c(FALSE, TRUE, FALSE),
    x3 = c(FALSE, FALSE, TRUE)
  )

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
})

test_that("loglik handles all components in candidate set", {
  # Create data where candidate set contains all components
  df <- data.frame(
    t = c(1, 2, 3),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, TRUE, TRUE),
    x2 = c(TRUE, TRUE, TRUE),
    x3 = c(TRUE, TRUE, TRUE)
  )

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
})

test_that("model works with custom column names", {
  # Create data with non-standard column names
  df <- data.frame(
    time = c(1, 2, 3),
    censor = c(TRUE, TRUE, TRUE),
    cand1 = c(TRUE, FALSE, FALSE),
    cand2 = c(FALSE, TRUE, TRUE),
    cand3 = c(TRUE, TRUE, TRUE)
  )

  model <- exp_series_md_c1_c2_c3(
    lifetime = "time",
    indicator = "censor",
    candset = "cand"
  )

  ll_fn <- loglik(model)
  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
})

test_that("model handles two-component systems", {
  df <- data.frame(
    t = c(1, 2, 3),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, TRUE, FALSE),
    x2 = c(FALSE, TRUE, TRUE)
  )

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)
  hess_fn <- hess_loglik(model)

  par <- c(0.5, 0.3)

  ll <- ll_fn(df, par)
  s <- score_fn(df, par)
  H <- hess_fn(df, par)

  expect_true(is.finite(ll))
  expect_length(s, 2)
  expect_equal(dim(H), c(2, 2))
})

test_that("model handles single observation", {
  df <- data.frame(
    t = 1.5,
    delta = TRUE,
    x1 = TRUE,
    x2 = FALSE,
    x3 = TRUE
  )

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
})


# ==============================================================================
# Tests for rdata method
# ==============================================================================

test_that("rdata.exp_series_md_c1_c2_c3 returns a function", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  expect_type(rdata_fn, "closure")
})

test_that("rdata generates data frame with correct structure", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  set.seed(42)
  df <- rdata_fn(theta = c(0.5, 0.3, 0.2), n = 100)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 100)
  expect_true("t" %in% names(df))
  expect_true("delta" %in% names(df))
  expect_true("x1" %in% names(df))
  expect_true("x2" %in% names(df))
  expect_true("x3" %in% names(df))
})

test_that("rdata respects custom column names from model", {
  model <- exp_series_md_c1_c2_c3(
    lifetime = "time",
    indicator = "censored",
    candset = "cand"
  )
  rdata_fn <- rdata(model)

  set.seed(42)
  df <- rdata_fn(theta = c(0.5, 0.3), n = 50)

  expect_true("time" %in% names(df))
  expect_true("censored" %in% names(df))
  expect_true("cand1" %in% names(df))
  expect_true("cand2" %in% names(df))
})

test_that("rdata applies right-censoring correctly", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  set.seed(42)
  df <- rdata_fn(theta = c(0.5, 0.3, 0.2), n = 200, tau = 1)

  # Check that censored observations exist and have correct properties
  censored <- !df$delta
  if (any(censored)) {
    # Censored times should equal tau
    expect_true(all(df$t[censored] == 1))
    # Censored observations should have empty candidate sets
    expect_true(all(!df$x1[censored]))
    expect_true(all(!df$x2[censored]))
    expect_true(all(!df$x3[censored]))
  }

  # Exact observations should have at least one component in candidate set
  exact <- df$delta
  if (any(exact)) {
    cand_sums <- df$x1[exact] + df$x2[exact] + df$x3[exact]
    expect_true(all(cand_sums >= 1))
  }
})

test_that("rdata generates candidate sets satisfying C1", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  # With p = 0, only the failed component should be in candidate set
  set.seed(42)
  df <- rdata_fn(theta = c(0.5, 0.3, 0.2), n = 100, p = 0)

  # Each exact observation should have exactly one component in candidate set
  exact <- df$delta
  cand_sums <- df$x1[exact] + df$x2[exact] + df$x3[exact]
  expect_true(all(cand_sums == 1))
})

test_that("rdata masking probability p works correctly", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  # With p = 1, all components should be in candidate set for exact obs
  set.seed(42)
  df <- rdata_fn(theta = c(0.5, 0.3, 0.2), n = 100, tau = Inf, p = 1)

  exact <- df$delta
  expect_true(all(df$x1[exact]))
  expect_true(all(df$x2[exact]))
  expect_true(all(df$x3[exact]))
})

test_that("rdata errors on invalid parameters", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  expect_error(rdata_fn(theta = c(-0.5, 0.3, 0.2), n = 100))
  expect_error(rdata_fn(theta = c(0, 0.3, 0.2), n = 100))
})

test_that("rdata generated data can be used with loglik", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)
  ll_fn <- loglik(model)

  set.seed(42)
  df <- rdata_fn(theta = c(0.5, 0.3, 0.2), n = 100, tau = 10, p = 0.3)
  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
})


# ==============================================================================
# Tests for observed_info method
# ==============================================================================

test_that("observed_info returns a function", {
  model <- exp_series_md_c1_c2_c3()
  obs_info_fn <- observed_info(model)

  expect_type(obs_info_fn, "closure")
})

test_that("observed_info returns negative of hessian", {
  model <- exp_series_md_c1_c2_c3()
  hess_fn <- hess_loglik(model)
  obs_info_fn <- observed_info(model)

  df <- create_exp_test_data(n = 50)
  par <- c(0.5, 0.3, 0.2)

  H <- hess_fn(df, par)
  I_obs <- obs_info_fn(df, par)

  expect_equal(I_obs, -H)
})

test_that("observed_info is positive semi-definite at MLE", {
  true_rates <- c(0.5, 0.3, 0.2)
  set.seed(42)
  df <- create_exp_test_data(n = 200, rates = true_rates)

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  obs_info_fn <- observed_info(model)

  # Find MLE
  result <- optim(
    par = c(1, 1, 1),
    fn = function(theta) -ll_fn(df, theta),
    method = "L-BFGS-B",
    lower = rep(1e-6, 3)
  )
  mle <- result$par

  I_obs <- obs_info_fn(df, par = mle)

  # Observed info should be positive semi-definite
  eigenvalues <- eigen(I_obs, symmetric = TRUE)$values
  expect_true(all(eigenvalues >= -1e-6))
})


# ==============================================================================
# Tests for fim method (Monte Carlo estimation)
# ==============================================================================

test_that("fim returns a function", {
  model <- exp_series_md_c1_c2_c3()
  fim_fn <- fim(model)

  expect_type(fim_fn, "closure")
})

test_that("fim returns matrix of correct dimensions", {
  model <- exp_series_md_c1_c2_c3()
  fim_fn <- fim(model)

  set.seed(42)
  I <- fim_fn(theta = c(0.5, 0.3, 0.2), n_obs = 50, n_samples = 50)

  expect_true(is.matrix(I))
  expect_equal(dim(I), c(3, 3))
  expect_true(all(is.finite(I)))
})

test_that("fim is approximately symmetric", {
  model <- exp_series_md_c1_c2_c3()
  fim_fn <- fim(model)

  set.seed(42)
  I <- fim_fn(theta = c(0.5, 0.3, 0.2), n_obs = 100, n_samples = 100)

  # MC FIM should be approximately symmetric
  expect_equal(I, t(I), tolerance = 0.1)
})

test_that("fim is approximately positive semi-definite", {
  model <- exp_series_md_c1_c2_c3()
  fim_fn <- fim(model)

  set.seed(42)
  I <- fim_fn(theta = c(0.5, 0.3, 0.2), n_obs = 100, n_samples = 200)

  # Symmetrize before eigenvalue check
  I_sym <- (I + t(I)) / 2
  eigenvalues <- eigen(I_sym, symmetric = TRUE)$values

  # All eigenvalues should be non-negative (allowing small numerical error)
  expect_true(all(eigenvalues >= -0.1))
})

test_that("fim accepts additional DGP parameters", {
  model <- exp_series_md_c1_c2_c3()
  fim_fn <- fim(model)

  set.seed(42)
  I <- fim_fn(theta = c(0.5, 0.3, 0.2), n_obs = 50, n_samples = 50,
              tau = 5, p = 0.3)

  expect_true(is.matrix(I))
  expect_equal(dim(I), c(3, 3))
})


# ==============================================================================
# Test: Cross-model consistency — Weibull(shape=1) ≡ Exponential
# ==============================================================================

test_that("Weibull(shape=1) loglik equals Exponential loglik", {
  # When all Weibull shapes = 1, Weibull(1, beta) = Exp(rate = 1/beta)
  set.seed(42)
  rates <- c(0.5, 0.3, 0.2)
  scales <- 1 / rates  # beta = 1/lambda

  # Generate data using exponential model
  exp_model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(exp_model)
  df <- rdata_fn(theta = rates, n = 100, tau = 10, p = 0.3)

  # Evaluate loglik with both models
  exp_ll <- loglik(exp_model)
  wei_model <- wei_series_md_c1_c2_c3()
  wei_ll <- loglik(wei_model)

  ll_exp <- exp_ll(df, par = rates)
  # Weibull par = (shape1, scale1, shape2, scale2, ...)
  wei_par <- as.numeric(rbind(rep(1, 3), scales))
  ll_wei <- wei_ll(df, par = wei_par)

  expect_equal(ll_exp, ll_wei, tolerance = 1e-8)
})

test_that("Weibull(shape=1) score is consistent with Exponential score", {
  set.seed(42)
  rates <- c(0.5, 0.3)
  scales <- 1 / rates

  exp_model <- exp_series_md_c1_c2_c3()
  df <- rdata(exp_model)(theta = rates, n = 50, tau = 10, p = 0.3)

  exp_score <- score(exp_model)(df, par = rates)
  wei_model <- wei_series_md_c1_c2_c3()
  wei_par <- c(1, scales[1], 1, scales[2])
  wei_score_val <- score(wei_model)(df, par = wei_par)

  # Weibull scale score relates to exp rate score:
  # d/d(beta) = d/d(lambda) * d(lambda)/d(beta) = d/d(lambda) * (-1/beta^2)
  # So exp score = -wei_scale_score * beta^2
  wei_scale_scores <- wei_score_val[c(2, 4)]
  exp_from_wei <- -wei_scale_scores * scales^2

  expect_equal(exp_score, exp_from_wei, tolerance = 1e-4)
})


# ==============================================================================
# Test: Single-component system (m=1)
# ==============================================================================

test_that("loglik works with single-component system (m=1)", {
  df <- data.frame(
    t = c(1, 2, 3),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, TRUE, TRUE)
  )

  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)
  hess_fn <- hess_loglik(model)

  par <- 0.5

  ll <- ll_fn(df, par = par)
  s <- score_fn(df, par = par)
  H <- hess_fn(df, par = par)

  expect_true(is.finite(ll))
  expect_length(s, 1)
  expect_equal(dim(H), c(1, 1))

  # For m=1 with no masking, loglik = -sum(t)*lambda + n*log(lambda)
  expected_ll <- -sum(df$t) * par + sum(df$delta) * log(par)
  expect_equal(ll, expected_ll, tolerance = 1e-10)
})

test_that("rdata works with single-component system (m=1)", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  set.seed(42)
  df <- rdata_fn(theta = 0.5, n = 50, tau = 10)

  expect_true("x1" %in% names(df))
  expect_true(all(df$x1[df$delta]))
})


# ==============================================================================
# Test: High censoring fraction
# ==============================================================================

test_that("loglik is finite with high censoring (>80%)", {
  model <- exp_series_md_c1_c2_c3()
  rdata_fn <- rdata(model)

  # Very low tau to get high censoring
  set.seed(42)
  df <- rdata_fn(theta = c(0.5, 0.3, 0.2), n = 200, tau = 0.3, p = 0.3)
  cens_frac <- mean(!df$delta)
  # Verify we actually have high censoring
  expect_gt(cens_frac, 0.5)

  ll_fn <- loglik(model)
  score_fn <- score(model)
  hess_fn <- hess_loglik(model)

  ll <- ll_fn(df, par = c(0.5, 0.3, 0.2))
  s <- score_fn(df, par = c(0.5, 0.3, 0.2))
  H <- hess_fn(df, par = c(0.5, 0.3, 0.2))

  expect_true(is.finite(ll))
  expect_true(all(is.finite(s)))
  expect_true(all(is.finite(H)))
})


# ==============================================================================
# Test: Error message validation
# ==============================================================================

test_that("error messages are informative", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)
  score_fn <- score(model)

  df_empty <- data.frame(t = numeric(0), x1 = logical(0))
  expect_error(ll_fn(df_empty, par = 0.5), "df is empty")

  df_no_t <- data.frame(x1 = c(TRUE, FALSE))
  expect_error(ll_fn(df_no_t, par = 0.5), "lifetime variable")

  df_no_cand <- data.frame(t = c(1, 2))
  expect_error(ll_fn(df_no_cand, par = 0.5), "no candidate set")

  # Wrong number of parameters (new standardized message)
  df_mismatch <- data.frame(t = c(1, 2), delta = c(TRUE, TRUE),
                             x1 = c(TRUE, FALSE), x2 = c(FALSE, TRUE))
  expect_error(ll_fn(df_mismatch, par = c(0.5)), "Expected 2 parameters")

  # Wrong number of parameters
  df <- data.frame(t = c(1, 2), delta = c(TRUE, TRUE),
                   x1 = c(TRUE, FALSE), x2 = c(FALSE, TRUE))
  expect_error(score_fn(df, par = c(0.5)), "parameters")
})


# ==============================================================================
# Test: C1 violation detection
# ==============================================================================

test_that("loglik detects C1 violation: exact observation with empty candidate set", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # Create data with C1 violation: delta=TRUE but all candidate indicators FALSE
  df_violation <- data.frame(
    t = c(1, 2, 3),
    delta = c(TRUE, TRUE, TRUE),
    x1 = c(TRUE, FALSE, FALSE),  # row 2 and 3 violate C1
    x2 = c(FALSE, FALSE, FALSE),
    x3 = c(FALSE, FALSE, FALSE)
  )

  expect_error(ll_fn(df_violation, par = c(0.5, 0.3, 0.2)), "C1 violated")
})

test_that("loglik allows empty candidate set for censored observations", {
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # Empty candidate sets are valid for censored observations (delta=FALSE)
  df_valid <- data.frame(
    t = c(1, 2, 3),
    delta = c(TRUE, FALSE, FALSE),
    x1 = c(TRUE, FALSE, FALSE),  # rows 2,3 are censored, empty set OK
    x2 = c(FALSE, FALSE, FALSE),
    x3 = c(FALSE, FALSE, FALSE)
  )

  ll <- ll_fn(df_valid, par = c(0.5, 0.3, 0.2))
  expect_true(is.finite(ll))
})
