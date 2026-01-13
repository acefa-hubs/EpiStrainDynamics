# Tests for diagnose_model() function
# Uses same fixtures as test-metrics.R to avoid duplication

# ==============================================================================
# SETUP: Load fixtures
# ==============================================================================

test_that("fixtures are available for diagnose_model tests", {
  skip_if_not(file.exists(test_path("fixtures/fit_rw_single.rds")),
              "Fixtures not available")
  expect_true(file.exists(test_path("fixtures/fit_rw_single.rds")))
})

# Load fixtures (reuse from test-metrics.R)
fit_rw_single <- readRDS(test_path("fixtures/fit_rw_single.rds"))
fit_ps_single <- readRDS(test_path("fixtures/fit_ps_single.rds"))
fit_rw_multi <- readRDS(test_path("fixtures/fit_rw_multi.rds"))
fit_ps_multi <- readRDS(test_path("fixtures/fit_ps_multi.rds"))

# ==============================================================================
# TESTS: INPUT VALIDATION
# ==============================================================================

test_that("diagnose_model() validates input class", {
  expect_error(
    diagnose_model(list(fit = "not a fit")),
    "must be an EpiStrainDynamics.fit object"
  )

  expect_error(
    diagnose_model("not a model"),
    "must be an EpiStrainDynamics.fit object"
  )

  expect_error(
    diagnose_model(data.frame(x = 1, y = 2)),
    "must be an EpiStrainDynamics.fit object"
  )
})

test_that("diagnose_model() validates threshold parameters", {
  # rhat_threshold should be numeric
  expect_error(
    diagnose_model(fit_rw_single, rhat_threshold = "bad"),
    "must be numeric"
  )

  # eff_sample_threshold should be numeric
  expect_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = "bad"),
    "must be numeric"
  )

  # rhat_threshold must be >= 1.0
  expect_error(
    diagnose_model(fit_rw_single, rhat_threshold = 0.5),
    "must be greater than or equal to 1.0"
  )

  expect_error(
    diagnose_model(fit_rw_single, rhat_threshold = 0.99),
    "must be greater than or equal to 1.0"
  )

  # eff_sample_threshold must be positive
  expect_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = 0),
    "must be positive"
  )

  expect_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = -100),
    "must be positive"
  )

  # Should accept valid values
  expect_no_error(
    diagnose_model(fit_rw_single, rhat_threshold = 1.0)
  )

  expect_no_error(
    diagnose_model(fit_rw_single, rhat_threshold = 1.1)
  )

  expect_no_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = 100)
  )
})

test_that("diagnose_model() warns about unusual threshold values", {
  # Very high rhat_threshold should warn
  expect_warning(
    diagnose_model(fit_rw_single, rhat_threshold = 5.0),
    "unusually high"
  )

  # Very low eff_sample_threshold should warn
  expect_warning(
    diagnose_model(fit_rw_single, eff_sample_threshold = 10),
    "very low"
  )

  # Very high eff_sample_threshold should warn
  expect_warning(
    diagnose_model(fit_rw_single, eff_sample_threshold = 10000),
    "very high"
  )
})

# ==============================================================================
# TESTS: OUTPUT STRUCTURE
# ==============================================================================

test_that("diagnose_model() returns correct structure", {
  diag <- diagnose_model(fit_rw_single)

  # Should return a list
  expect_type(diag, "list")

  # Should have expected components
  expect_named(diag, c("convergence", "rhat_issues", "eff_sample_issues",
                       "max_rhat", "min_neff", "summary"))
})

test_that("diagnose_model() convergence is logical", {
  diag <- diagnose_model(fit_rw_single)

  expect_type(diag$convergence, "logical")
  expect_length(diag$convergence, 1)
  expect_false(is.na(diag$convergence))
})

test_that("diagnose_model() issues are character vectors", {
  diag <- diagnose_model(fit_rw_single)

  expect_type(diag$rhat_issues, "character")
  expect_type(diag$eff_sample_issues, "character")

  # Can be empty (no issues)
  expect_true(length(diag$rhat_issues) >= 0)
  expect_true(length(diag$eff_sample_issues) >= 0)
})

test_that("diagnose_model() max_rhat and min_neff are numeric", {
  diag <- diagnose_model(fit_rw_single)

  expect_type(diag$max_rhat, "double")
  expect_type(diag$min_neff, "double")

  expect_length(diag$max_rhat, 1)
  expect_length(diag$min_neff, 1)

  expect_true(is.finite(diag$max_rhat))
  expect_true(is.finite(diag$min_neff))
})

test_that("diagnose_model() summary is a matrix", {
  diag <- diagnose_model(fit_rw_single)

  expect_true(is.matrix(diag$summary))

  # Should have standard Stan summary columns
  expected_cols <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%",
                     "75%", "97.5%", "n_eff", "Rhat")
  expect_true(all(expected_cols %in% colnames(diag$summary)))
})

# ==============================================================================
# TESTS: CONVERGENCE LOGIC
# ==============================================================================

test_that("diagnose_model() convergence logic is correct", {
  diag <- diagnose_model(fit_rw_single, rhat_threshold = 1.1,
                         eff_sample_threshold = 100)

  # If convergence is TRUE, there should be no issues
  if (diag$convergence) {
    expect_equal(length(diag$rhat_issues), 0)
    expect_equal(length(diag$eff_sample_issues), 0)
  }

  # If convergence is FALSE, there should be at least one issue
  if (!diag$convergence) {
    total_issues <- length(diag$rhat_issues) + length(diag$eff_sample_issues)
    expect_gt(total_issues, 0)
  }
})

# ==============================================================================
# TESTS: WORKS WITH DIFFERENT MODEL TYPES
# ==============================================================================

test_that("diagnose_model() works with single pathogen random walk", {
  expect_no_error(diag <- diagnose_model(fit_rw_single))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

test_that("diagnose_model() works with single pathogen p-spline", {
  expect_no_error(diag <- diagnose_model(fit_ps_single))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

test_that("diagnose_model() works with multiple pathogen random walk", {
  expect_no_error(diag <- diagnose_model(fit_rw_multi))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

test_that("diagnose_model() works with multiple pathogen p-spline", {
  expect_no_error(diag <- diagnose_model(fit_ps_multi))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

# ==============================================================================
# TESTS: DIAGNOSTIC VALUES ARE REASONABLE
# ==============================================================================

test_that("diagnose_model() max_rhat is >= 1", {
  diag <- diagnose_model(fit_rw_single)

  # R-hat should always be >= 1 by definition
  expect_gte(diag$max_rhat, 1.0)
})

test_that("diagnose_model() min_neff is positive", {
  diag <- diagnose_model(fit_rw_single)

  # Effective sample size should be positive
  expect_gt(diag$min_neff, 0)
})

test_that("diagnose_model() identifies parameter names in issues", {
  diag <- diagnose_model(fit_rw_single, rhat_threshold = 1.01,
                         eff_sample_threshold = 500)

  # If there are issues, they should be named parameters
  if (length(diag$rhat_issues) > 0) {
    expect_true(all(nzchar(diag$rhat_issues)))

    # Should be valid row names from summary
    expect_true(all(diag$rhat_issues %in% rownames(diag$summary)))
  }

  if (length(diag$eff_sample_issues) > 0) {
    expect_true(all(nzchar(diag$eff_sample_issues)))
    expect_true(all(diag$eff_sample_issues %in% rownames(diag$summary)))
  }
})

# ==============================================================================
# TESTS: CONSOLE OUTPUT (captured)
# ==============================================================================

test_that("diagnose_model() prints summary to console", {
  output <- capture.output(diag <- diagnose_model(fit_rw_single))

  # Should have printed something
  expect_gt(length(output), 0)

  # Should contain key phrases
  expect_true(any(grepl("Model Convergence Diagnostics", output)))
  expect_true(any(grepl("Overall convergence", output)))
  expect_true(any(grepl("Maximum R-hat", output)))
  expect_true(any(grepl("Minimum n_eff", output)))
})

test_that("diagnose_model() shows GOOD or ISSUES DETECTED", {
  output <- capture.output(diag <- diagnose_model(fit_rw_single))

  # Should contain either GOOD or ISSUES DETECTED
  output_text <- paste(output, collapse = " ")
  expect_true(grepl("GOOD|ISSUES DETECTED", output_text))
})

test_that("diagnose_model() console output matches convergence status", {
  output <- capture.output(diag <- diagnose_model(fit_rw_single))
  output_text <- paste(output, collapse = " ")

  if (diag$convergence) {
    expect_true(grepl("GOOD", output_text))
  } else {
    expect_true(grepl("ISSUES DETECTED", output_text))
  }
})

# ==============================================================================
# TESTS: INVISIBLE RETURN
# ==============================================================================

test_that("diagnose_model() returns invisibly", {
  # The function should return invisible
  # We can't directly test invisibility, but we can check that
  # assigning doesn't print
  output <- capture.output({
    result <- diagnose_model(fit_rw_single)
  })

  # Should have printed diagnostics but not returned value
  expect_true(any(grepl("Model Convergence", output)))
})

# ==============================================================================
# TESTS: COMPARISON ACROSS MODEL TYPES
# ==============================================================================

test_that("diagnose_model() produces consistent output across model types", {
  diag_rw_single <- diagnose_model(fit_rw_single)
  diag_ps_single <- diagnose_model(fit_ps_single)
  diag_rw_multi <- diagnose_model(fit_rw_multi)
  diag_ps_multi <- diagnose_model(fit_ps_multi)

  # All should have the same structure
  expected_names <- c("convergence", "rhat_issues", "eff_sample_issues",
                      "max_rhat", "min_neff", "summary")

  expect_named(diag_rw_single, expected_names)
  expect_named(diag_ps_single, expected_names)
  expect_named(diag_rw_multi, expected_names)
  expect_named(diag_ps_multi, expected_names)
})

test_that("diagnose_model() summary has more parameters for multi-pathogen models", {
  diag_single <- diagnose_model(fit_rw_single)
  diag_multi <- diagnose_model(fit_rw_multi)

  # Multi-pathogen model should have more parameters
  n_params_single <- nrow(diag_single$summary)
  n_params_multi <- nrow(diag_multi$summary)

  expect_gt(n_params_multi, n_params_single)
})

# ==============================================================================
# TESTS: COMPARISON WITH RSTAN DIAGNOSTICS
# ==============================================================================

test_that("diagnose_model() max_rhat matches Stan summary", {
  diag <- diagnose_model(fit_rw_single)

  # Get R-hat values directly from Stan
  stan_summary <- rstan::summary(fit_rw_single$fit)$summary
  stan_rhat <- stan_summary[, "Rhat"]
  stan_max_rhat <- max(stan_rhat, na.rm = TRUE)

  # Should match
  expect_equal(diag$max_rhat, stan_max_rhat)
})

test_that("diagnose_model() min_neff matches Stan summary", {
  diag <- diagnose_model(fit_rw_single)

  # Get n_eff values directly from Stan
  stan_summary <- rstan::summary(fit_rw_single$fit)$summary
  stan_neff <- stan_summary[, "n_eff"]
  stan_min_neff <- min(stan_neff, na.rm = TRUE)

  # Should match
  expect_equal(diag$min_neff, stan_min_neff)
})

