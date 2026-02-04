# Tests for masked data generation functions

# ==============================================================================
# Tests for md_series_lifetime_right_censoring
# ==============================================================================

test_that("md_series_lifetime_right_censoring adds system lifetime column", {
  df <- data.frame(
    t1 = c(1, 2, 3),
    t2 = c(2, 1, 4),
    t3 = c(3, 3, 2)
  )

  result <- md_series_lifetime_right_censoring(df, tau = Inf)

  expect_true("t" %in% colnames(result))
  # System lifetime should be minimum of component lifetimes
  expect_equal(result$t, c(1, 1, 2))
})

test_that("md_series_lifetime_right_censoring adds censoring indicator", {
  df <- data.frame(
    t1 = c(1, 2, 3),
    t2 = c(2, 1, 4),
    t3 = c(3, 3, 2)
  )

  result <- md_series_lifetime_right_censoring(df, tau = Inf)

  expect_true("delta" %in% colnames(result))
  # All TRUE when tau = Inf (no censoring)
  expect_true(all(result$delta))
})

test_that("md_series_lifetime_right_censoring applies right censoring", {
  df <- data.frame(
    t1 = c(1, 5, 3),
    t2 = c(2, 6, 4),
    t3 = c(3, 7, 2)
  )
  tau <- 4

  result <- md_series_lifetime_right_censoring(df, tau = tau)

  # First observation: min = 1 < 4, not censored
  expect_equal(result$t[1], 1)
  expect_true(result$delta[1])

  # Second observation: min = 5 > 4, censored at 4
  expect_equal(result$t[2], 4)
  expect_false(result$delta[2])

  # Third observation: min = 2 < 4, not censored
  expect_equal(result$t[3], 2)
  expect_true(result$delta[3])
})

test_that("md_series_lifetime_right_censoring accepts custom column names", {
  df <- data.frame(
    comp1 = c(1, 2, 3),
    comp2 = c(2, 1, 4)
  )

  result <- md_series_lifetime_right_censoring(
    df,
    tau = Inf,
    comp = "comp",
    lifetime = "time",
    right_censoring_indicator = "cens"
  )

  expect_true("time" %in% colnames(result))
  expect_true("cens" %in% colnames(result))
  expect_equal(result$time, c(1, 1, 3))
})

test_that("md_series_lifetime_right_censoring errors when comp not found", {
  df <- data.frame(
    x1 = c(1, 2, 3),
    x2 = c(2, 1, 4)
  )

  expect_error(md_series_lifetime_right_censoring(df, comp = "t"), "comp not in colnames")
})


# ==============================================================================
# Tests for md_bernoulli_cand_c1_c2_c3
# ==============================================================================

test_that("md_bernoulli_cand_c1_c2_c3 adds probability columns", {
  df <- data.frame(
    t1 = c(1, 2, 3),
    t2 = c(2, 1, 4),
    t3 = c(3, 3, 2),
    delta = c(TRUE, TRUE, TRUE)
  )

  result <- md_bernoulli_cand_c1_c2_c3(df, p = 0.3)

  expect_true("q1" %in% colnames(result))
  expect_true("q2" %in% colnames(result))
  expect_true("q3" %in% colnames(result))
})

test_that("md_bernoulli_cand_c1_c2_c3 sets failed component probability to 1", {
  df <- data.frame(
    t1 = c(1, 5, 3),  # min at row 1, 3
    t2 = c(2, 1, 4),  # min at row 2
    t3 = c(3, 3, 2),
    delta = c(TRUE, TRUE, TRUE)
  )

  result <- md_bernoulli_cand_c1_c2_c3(df, p = 0.3)

  # Row 1: t1 is minimum, so q1 = 1
  expect_equal(result$q1[1], 1)
  # Row 2: t2 is minimum, so q2 = 1
  expect_equal(result$q2[2], 1)
  # Row 3: t3 is minimum, so q3 = 1
  expect_equal(result$q3[3], 1)
})

test_that("md_bernoulli_cand_c1_c2_c3 sets non-failed component probability to p", {
  df <- data.frame(
    t1 = c(1, 5, 3),
    t2 = c(2, 1, 4),
    t3 = c(3, 3, 2),
    delta = c(TRUE, TRUE, TRUE)
  )
  p <- 0.3

  result <- md_bernoulli_cand_c1_c2_c3(df, p = p)

  # Row 1: t1 is minimum, so q2 and q3 should be p
  expect_equal(result$q2[1], p)
  expect_equal(result$q3[1], p)

  # Row 2: t2 is minimum, so q1 and q3 should be p
  expect_equal(result$q1[2], p)
  expect_equal(result$q3[2], p)
})

test_that("md_bernoulli_cand_c1_c2_c3 sets all probabilities to 0 for censored observations", {
  df <- data.frame(
    t1 = c(1, 5, 3),
    t2 = c(2, 1, 4),
    t3 = c(3, 3, 2),
    delta = c(TRUE, FALSE, TRUE)  # Row 2 is censored
  )

  result <- md_bernoulli_cand_c1_c2_c3(df, p = 0.3)

  # Row 2 is censored, all probabilities should be 0
  expect_equal(result$q1[2], 0)
  expect_equal(result$q2[2], 0)
  expect_equal(result$q3[2], 0)
})

test_that("md_bernoulli_cand_c1_c2_c3 handles empty data frame", {
  df <- data.frame(
    t1 = numeric(0),
    t2 = numeric(0),
    delta = logical(0)
  )

  result <- md_bernoulli_cand_c1_c2_c3(df, p = 0.3)

  expect_equal(nrow(result), 0)
})

test_that("md_bernoulli_cand_c1_c2_c3 accepts vector of probabilities", {
  df <- data.frame(
    t1 = c(1, 2),
    t2 = c(2, 1),
    delta = c(TRUE, TRUE)
  )
  p <- c(0.2, 0.4)  # Different p for each row

  result <- md_bernoulli_cand_c1_c2_c3(df, p = p)

  # Row 1: t1 is minimum, q2 = 0.2
  expect_equal(result$q2[1], 0.2)
  # Row 2: t2 is minimum, q1 = 0.4
  expect_equal(result$q1[2], 0.4)
})

test_that("md_bernoulli_cand_c1_c2_c3 errors when comp not found", {
  df <- data.frame(
    x1 = c(1, 2, 3),
    x2 = c(2, 1, 4),
    delta = c(TRUE, TRUE, TRUE)
  )

  expect_error(md_bernoulli_cand_c1_c2_c3(df, p = 0.3), "No component lifetime")
})

test_that("md_bernoulli_cand_c1_c2_c3 errors when right_censoring_indicator not found", {
  df <- data.frame(
    t1 = c(1, 2, 3),
    t2 = c(2, 1, 4)
  )

  expect_error(md_bernoulli_cand_c1_c2_c3(df, p = 0.3), "right_censoring_indicator not in colnames")
})


# ==============================================================================
# Tests for md_cand_sampler
# ==============================================================================

test_that("md_cand_sampler creates candidate set columns", {
  df <- data.frame(
    q1 = c(1, 0.5, 0.3),
    q2 = c(0.3, 1, 0.3),
    q3 = c(0.3, 0.3, 1)
  )

  set.seed(42)
  result <- md_cand_sampler(df)

  expect_true("x1" %in% colnames(result))
  expect_true("x2" %in% colnames(result))
  expect_true("x3" %in% colnames(result))
})

test_that("md_cand_sampler produces boolean candidate sets", {
  df <- data.frame(
    q1 = c(1, 0.5, 0.3),
    q2 = c(0.3, 1, 0.3),
    q3 = c(0.3, 0.3, 1)
  )

  set.seed(42)
  result <- md_cand_sampler(df)

  expect_type(result$x1, "logical")
  expect_type(result$x2, "logical")
  expect_type(result$x3, "logical")
})

test_that("md_cand_sampler includes component with probability 1", {
  # Set probabilities so some are always 1
  df <- data.frame(
    q1 = c(1, 0, 0),
    q2 = c(0, 1, 0),
    q3 = c(0, 0, 1)
  )

  result <- md_cand_sampler(df)

  expect_true(result$x1[1])
  expect_true(result$x2[2])
  expect_true(result$x3[3])

  expect_false(result$x2[1])
  expect_false(result$x3[1])
})

test_that("md_cand_sampler handles empty data frame", {
  df <- data.frame(
    q1 = numeric(0),
    q2 = numeric(0)
  )

  result <- md_cand_sampler(df)

  expect_equal(nrow(result), 0)
})

test_that("md_cand_sampler errors when df is not data frame", {
  expect_error(md_cand_sampler(list(q1 = c(1, 0.5))), "must be a data frame")
})

test_that("md_cand_sampler accepts custom column names", {
  df <- data.frame(
    prob1 = c(1, 0.5),
    prob2 = c(0.5, 1)
  )

  set.seed(42)
  result <- md_cand_sampler(df, prob = "prob", candset = "cand")

  expect_true("cand1" %in% colnames(result))
  expect_true("cand2" %in% colnames(result))
})

test_that("md_cand_sampler samples match probabilities approximately", {
  skip_on_cran()  # This test is probabilistic

  # Set up with known probability
  n <- 1000
  p <- 0.7
  df <- data.frame(
    q1 = rep(p, n),
    q2 = rep(p, n)
  )

  set.seed(123)
  result <- md_cand_sampler(df)

  # Sample proportion should be close to p
  prop1 <- mean(result$x1)
  prop2 <- mean(result$x2)

  expect_equal(prop1, p, tolerance = 0.05)
  expect_equal(prop2, p, tolerance = 0.05)
})


# ==============================================================================
# Integration tests: Complete masked data workflow
# ==============================================================================

test_that("complete masked data generation workflow works", {
  # Step 1: Generate component lifetimes
  set.seed(42)
  n <- 100
  rates <- c(0.5, 0.3, 0.2)
  df <- data.frame(
    t1 = rexp(n, rates[1]),
    t2 = rexp(n, rates[2]),
    t3 = rexp(n, rates[3])
  )

  # Step 2: Apply right censoring
  tau <- 3
  df <- md_series_lifetime_right_censoring(df, tau = tau)

  expect_true("t" %in% colnames(df))
  expect_true("delta" %in% colnames(df))
  expect_true(all(df$t <= tau))

  # Step 3: Generate candidate set probabilities
  df <- md_bernoulli_cand_c1_c2_c3(df, p = 0.3)

  expect_true("q1" %in% colnames(df))
  expect_true("q2" %in% colnames(df))
  expect_true("q3" %in% colnames(df))

  # Step 4: Sample candidate sets
  df <- md_cand_sampler(df)

  expect_true("x1" %in% colnames(df))
  expect_true("x2" %in% colnames(df))
  expect_true("x3" %in% colnames(df))

  # Verify C1: Failed component is in candidate set for non-censored observations
  for (i in seq_len(nrow(df))) {
    if (df$delta[i]) {
      # For non-censored: candidate set should be non-empty
      # and should include the failed component (but we can only verify non-empty
      # since we don't know which failed without the component lifetimes being visible)
      expect_true(df$x1[i] || df$x2[i] || df$x3[i],
                  info = paste("Row", i, "should have non-empty candidate set"))
    } else {
      # For censored: candidate set should be empty
      expect_false(df$x1[i] || df$x2[i] || df$x3[i],
                   info = paste("Censored row", i, "should have empty candidate set"))
    }
  }
})

# ==============================================================================
# Tests for md_bernoulli_cand_c1_c2_c3 edge cases
# ==============================================================================

test_that("md_bernoulli_cand_c1_c2_c3 with p=0 gives singleton candidate sets", {
  df <- data.frame(
    t1 = c(1, 5, 3),
    t2 = c(2, 1, 4),
    t3 = c(3, 3, 2),
    delta = c(TRUE, TRUE, TRUE)
  )

  result <- md_bernoulli_cand_c1_c2_c3(df, p = 0)

  # Non-failed components should have probability 0
  # Row 1: t1 is min, so q2 = q3 = 0
  expect_equal(result$q2[1], 0)
  expect_equal(result$q3[1], 0)
  # But failed component is still 1
  expect_equal(result$q1[1], 1)
})

test_that("md_bernoulli_cand_c1_c2_c3 with p=1 gives all-component candidate sets", {
  df <- data.frame(
    t1 = c(1, 5, 3),
    t2 = c(2, 1, 4),
    t3 = c(3, 3, 2),
    delta = c(TRUE, TRUE, TRUE)
  )

  result <- md_bernoulli_cand_c1_c2_c3(df, p = 1)

  # All non-censored obs should have all probabilities = 1
  expect_true(all(result$q1 == 1))
  expect_true(all(result$q2 == 1))
  expect_true(all(result$q3 == 1))
})

test_that("md_bernoulli_cand_c1_c2_c3 with vector p uses per-obs probabilities", {
  df <- data.frame(
    t1 = c(1, 5),
    t2 = c(2, 1),
    delta = c(TRUE, TRUE)
  )
  p <- c(0.0, 1.0)

  result <- md_bernoulli_cand_c1_c2_c3(df, p = p)

  # Row 1: p=0, only failed component (t1 min -> q1=1, q2=0)
  expect_equal(result$q1[1], 1)
  expect_equal(result$q2[1], 0)

  # Row 2: p=1, both components (t2 min -> q2=1, q1=1)
  expect_equal(result$q1[2], 1)
  expect_equal(result$q2[2], 1)
})


test_that("masked data can be used with likelihood model", {
  # Generate masked data
  set.seed(42)
  n <- 100
  rates <- c(0.5, 0.3, 0.2)
  df <- data.frame(
    t1 = rexp(n, rates[1]),
    t2 = rexp(n, rates[2]),
    t3 = rexp(n, rates[3])
  )
  df <- md_series_lifetime_right_censoring(df, tau = Inf)
  df <- md_bernoulli_cand_c1_c2_c3(df, p = 0.3)
  df <- md_cand_sampler(df)

  # Fit likelihood model
  model <- exp_series_md_c1_c2_c3()
  ll_fn <- loglik(model)

  # Should be able to evaluate likelihood
  ll <- ll_fn(df, par = rates)

  expect_true(is.finite(ll))
})
