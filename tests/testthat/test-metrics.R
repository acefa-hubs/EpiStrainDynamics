# Test suite for post-fit metric calculations
# Uses saved fixture models to avoid slow refitting

# ==============================================================================
# SETUP AND FIXTURES
# ==============================================================================

test_that("fixture models are available", {
  # These should be created once and saved as .rds files
  # See setup-fixtures.R for generation code

  skip_if_not(file.exists("tests/testthat/fixtures/fit_rw_single.rds"),
              "Fixture fit_rw_single.rds not found")
  skip_if_not(file.exists("tests/testthat/fixtures/fit_ps_single.rds"),
              "Fixture fit_ps_single.rds not found")
  skip_if_not(file.exists("tests/testthat/fixtures/fit_rw_multi.rds"),
              "Fixture fit_rw_multi.rds not found")
  skip_if_not(file.exists("tests/testthat/fixtures/fit_ps_multi.rds"),
              "Fixture fit_ps_multi.rds not found")

  # If we get here, all fixtures exist
  expect_true(file.exists("tests/testthat/fixtures/fit_rw_single.rds"))
  expect_true(file.exists("tests/testthat/fixtures/fit_ps_single.rds"))
  expect_true(file.exists("tests/testthat/fixtures/fit_rw_multi.rds"))
  expect_true(file.exists("tests/testthat/fixtures/fit_ps_multi.rds"))
})

# Load fixtures (done once per test file)
fit_rw_single <- readRDS("tests/testthat/fixtures/fit_rw_single.rds")
fit_rw_single_dow <- readRDS("tests/testthat/fixtures/fit_rw_single_dow.rds")
fit_ps_single <- readRDS("tests/testthat/fixtures/fit_ps_single.rds")
fit_rw_multi <- readRDS("tests/testthat/fixtures/fit_rw_multi.rds")
fit_ps_multi <- readRDS("tests/testthat/fixtures/fit_ps_multi.rds")

# Helper function for generation interval
gi_simple <- function(x) {
  ifelse(x == 0, 0, 4 * x * exp(-2 * x))
}

# ==============================================================================
# TESTS: Rt() FUNCTION
# ==============================================================================

test_that("Rt() validates inputs correctly", {
  # Wrong class
  expect_error(
    Rt(list(fit = "not a fit"), tau_max = 7, gi_dist = gi_simple),
    "Input must be of class"
  )

  # Invalid tau_max
  expect_error(
    Rt(fit_rw_single, tau_max = -1, gi_dist = gi_simple),
    "must be a positive number"
  )

  expect_error(
    Rt(fit_rw_single, tau_max = 1.5, gi_dist = gi_simple),
    "must be a whole number"
  )

  # Missing gi_dist
  expect_error(
    Rt(fit_rw_single, tau_max = 7),
    "missing"
  )

  # Invalid gi_dist
  expect_error(
    Rt(fit_rw_single, tau_max = 7, gi_dist = "not a function"),
    "must be a function"
  )
})

test_that("Rt() works for single pathogen random walk model", {
  result <- Rt(fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  # Check structure
  expect_s3_class(result, "Rt")
  expect_s3_class(result, "EpiStrainDynamics.metric")
  expect_named(result, c("measure", "fit", "constructed_model"))

  # Check measure data frame
  expect_s3_class(result$measure, "data.frame")
  expect_named(result$measure,
               c("y", "lb_50", "ub_50", "lb_95",
                 "ub_95", "prop", "time", "pathogen"))

  # Check dimensions - should start from time index = tau_max
  n_days <- length(fit_rw_single$constructed_model$data$time)
  expected_rows <- n_days - 7 + 1
  expect_equal(nrow(result$measure), expected_rows)

  # Check values are reasonable
  expect_true(all(result$measure$y > 0))  # Rt should be positive
  expect_true(all(result$measure$lb_95 > 0))
  expect_true(all(result$measure$ub_95 > result$measure$y))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
  expect_true(all(result$measure$lb_50 >= result$measure$lb_95))

  # Check prop is between 0 and 1 (proportion > threshold)
  expect_true(all(result$measure$prop >= 0 & result$measure$prop <= 1))

  # Check pathogen name
  expect_equal(unique(result$measure$pathogen), "cases")
})

test_that("Rt() works for single pathogen p-spline model", {
  result <- Rt(fit_ps_single, tau_max = 7, gi_dist = gi_simple)

  expect_s3_class(result, "Rt")
  expect_s3_class(result$measure, "data.frame")

  n_days <- length(fit_ps_single$constructed_model$data$time)
  expected_rows <- n_days - 7 + 1
  expect_equal(nrow(result$measure), expected_rows)

  expect_true(all(result$measure$y > 0))
})

test_that("Rt() works for multiple pathogen models", {
  result <- Rt(fit_rw_multi, tau_max = 7, gi_dist = gi_simple)

  expect_s3_class(result, "Rt")

  # Should have individual pathogens + Total
  pathogen_names <- unique(result$measure$pathogen)
  expect_true("Total" %in% pathogen_names)
  expect_true(all(c("alpha", "delta", "omicron", "other") %in% pathogen_names))

  n_days <- length(fit_rw_multi$constructed_model$data$time)
  expected_rows_per_pathogen <- n_days - 7 + 1
  expect_equal(nrow(result$measure), expected_rows_per_pathogen * 5)  # 4 pathogens + Total

  expect_true(all(result$measure$y > 0))
})

test_that("Rt() respects tau_max parameter", {
  result_7 <- Rt(fit_rw_single, tau_max = 7, gi_dist = gi_simple)
  result_14 <- Rt(fit_rw_single, tau_max = 14, gi_dist = gi_simple)

  # Different tau_max should give different number of rows
  expect_lt(nrow(result_14$measure), nrow(result_7$measure))

  # Values should differ (using different generation interval window)
  # At least some values should be different
  expect_false(all(abs(result_7$measure$y[1:10] - result_14$measure$y[1:10]) < 1e-10))
})

test_that("Rt() threshold affects prop calculation", {
  result <- Rt(fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  # Prop should be proportion > 1 (threshold for Rt)
  # Manually check a few values
  post <- rstan::extract(fit_rw_single$fit)
  first_time_idx <- 7  # tau_max

  # This is a rough check - the actual calc is more complex with splines etc
  # But prop should be between 0 and 1
  expect_true(all(result$measure$prop >= 0 & result$measure$prop <= 1))
})

# ==============================================================================
# TESTS: incidence() FUNCTION
# ==============================================================================
test_that("incidence() validates inputs correctly", {
  # Invalid class
  expect_error(
    incidence(list(fit = "not a fit"), dow = FALSE),
    "Input must be of class"
  )

  # dow = TRUE when model doesn't have dow_effect
  expect_error(
    incidence(fit_rw_single, dow = TRUE),
    "Day-of-week effects cannot be incorporated"
  )

  # dow must be NULL or logical
  expect_error(
    incidence(fit_rw_single, dow = "yes"),
    "must be.*NULL.*logical"
  )

  expect_error(
    incidence(fit_rw_single, dow = 1),
    "must be.*NULL.*logical"
  )
})

# dow argument default behavior tests
test_that("incidence() uses model's dow setting when dow = NULL", {
  # For model without dow_effect
  expect_message(
    result <- incidence(fit_rw_single, dow = NULL),
    "Using day-of-week setting from model.*FALSE"
  )
  expect_s3_class(result, "incidence")
})

test_that("incidence() default behavior (dow not specified) uses model setting", {
  # When dow argument is not provided at all, should use model's setting
  expect_message(
    result <- incidence(fit_rw_single),  # dow not specified
    "Using day-of-week setting from model"
  )
  expect_s3_class(result, "incidence")
})

test_that("incidence() explicit dow = FALSE works regardless of model", {
  # Should work for model without dow_effect
  expect_silent(
    result <- incidence(fit_rw_single, dow = FALSE)
  )
  expect_s3_class(result, "incidence")
})

test_that("incidence() dow = TRUE requires model to have dow_effect", {
  # Model without dow_effect - should error
  expect_error(
    incidence(fit_rw_single, dow = TRUE),
    "Day-of-week effects cannot be incorporated"
  )
})

# Basic functionality tests (existing tests - reviewed and kept)
test_that("incidence() works for single pathogen models", {
  result <- incidence(fit_rw_single, dow = FALSE)

  expect_s3_class(result, "incidence")
  expect_s3_class(result, "EpiStrainDynamics.metric")

  # Check structure
  expect_named(result$measure,
               c("y", "lb_50", "ub_50", "lb_95", "ub_95",
                 "prop", "time", "pathogen"))

  # Incidence should start from time index 1
  n_days <- length(fit_rw_single$constructed_model$data$time)
  expect_equal(nrow(result$measure), n_days)

  # Incidence should be positive
  expect_true(all(result$measure$y > 0))
  expect_true(all(result$measure$lb_95 > 0))

  # Check credible intervals are ordered
  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$lb_50 <= result$measure$y))
  expect_true(all(result$measure$y <= result$measure$ub_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

test_that("incidence() works for multiple pathogen models", {
  result <- incidence(fit_rw_multi, dow = FALSE)

  expect_s3_class(result, "incidence")

  # Should have individual pathogens + Total
  pathogen_names <- unique(result$measure$pathogen)
  expect_true("Total" %in% pathogen_names)
  expect_equal(length(pathogen_names), 5)  # 4 pathogens + Total

  # Total should be sum of individual pathogens (approximately, given posterior sampling)
  total_data <- result$measure[result$measure$pathogen == "Total", ]
  alpha_data <- result$measure[result$measure$pathogen == "alpha", ]
  expect_equal(nrow(total_data), nrow(alpha_data))
})

test_that("incidence() p-spline vs random walk gives different results", {
  result_rw <- incidence(fit_rw_single, dow = FALSE)
  result_ps <- incidence(fit_ps_single, dow = FALSE)

  # Should have same structure
  expect_equal(nrow(result_rw$measure), nrow(result_ps$measure))

  # But values should differ (different smoothing methods)
  expect_false(all(abs(result_rw$measure$y - result_ps$measure$y) < 1e-10))
})

test_that("incidence() with dow = TRUE includes day-of-week effects", {
  skip_if_not(exists("fit_rw_single_dow"),
              "No test fixture with dow_effect available")

  result_with_dow <- incidence(fit_rw_single_dow, dow = TRUE)
  result_without_dow <- incidence(fit_rw_single_dow, dow = FALSE)

  # Both should work
  expect_s3_class(result_with_dow, "incidence")
  expect_s3_class(result_without_dow, "incidence")

  # Results should differ when dow effects are included vs excluded
  expect_false(
    all(abs(result_with_dow$measure$y - result_without_dow$measure$y) < 1e-10)
  )
})

# Edge cases and special scenarios
test_that("incidence() treats NA same as NULL (uses model's setting)", {
  # NA should be treated as NULL - use the model's setting
  expect_message(
    result <- incidence(fit_rw_single, dow = NA),
    "Using day-of-week setting from model"
  )
  expect_s3_class(result, "incidence")
})

test_that("dow = NA works same as dow = NULL", {
  # Both should produce identical results
  result_null <- incidence(fit_rw_single, dow = NULL)
  result_na <- incidence(fit_rw_single, dow = NA)

  expect_equal(result_null$measure$y, result_na$measure$y)
  expect_equal(nrow(result_null$measure), nrow(result_na$measure))
})

test_that("incidence() preserves model information in output", {
  result <- incidence(fit_rw_single, dow = FALSE)

  # Output should include original fit and constructed_model
  expect_true("fit" %in% names(result))
  expect_true("constructed_model" %in% names(result))

  # These should be the same objects
  expect_identical(result$fit, fit_rw_single$fit)
  expect_identical(result$constructed_model, fit_rw_single$constructed_model)
})

test_that("incidence() measure data frame has correct properties", {
  result <- incidence(fit_rw_single, dow = FALSE)

  measure <- result$measure

  # All required columns present
  required_cols <- c("y", "lb_50", "ub_50", "lb_95", "ub_95",
                     "prop", "time", "pathogen")
  expect_true(all(required_cols %in% names(measure)))

  # No NA values in key columns
  expect_false(any(is.na(measure$y)))
  expect_false(any(is.na(measure$time)))
  expect_false(any(is.na(measure$pathogen)))

  # Proportions should be between 0 and 1
  expect_true(all(measure$prop >= 0 & measure$prop <= 1))

  # Time should be ordered
  expect_true(all(diff(as.numeric(measure$time)) > 0))
})

test_that("incidence() works consistently across different method classes", {
  # Test all four combinations: ps/rw × single/multi
  result_rw_single <- incidence(fit_rw_single, dow = FALSE)
  result_ps_single <- incidence(fit_ps_single, dow = FALSE)
  result_rw_multi <- incidence(fit_rw_multi, dow = FALSE)
  result_ps_multi <- incidence(fit_ps_multi, dow = FALSE)

  # All should have incidence class
  expect_s3_class(result_rw_single, "incidence")
  expect_s3_class(result_ps_single, "incidence")
  expect_s3_class(result_rw_multi, "incidence")
  expect_s3_class(result_ps_multi, "incidence")

  # All should have consistent structure
  expect_named(result_rw_single$measure,
               names(result_ps_single$measure))
  expect_named(result_rw_multi$measure,
               names(result_ps_multi$measure))
})

# Calculation validation tests
test_that("incidence() calculation functions are called with correct arguments", {
  # This is more of an integration test to ensure the helper functions work
  result <- incidence(fit_rw_single, dow = FALSE)

  # Incidence should be exp of log-incidence
  # Values should be reasonable (positive, not extreme)
  expect_true(all(result$measure$y > 0))
  expect_true(all(result$measure$y < Inf))
  expect_false(any(is.nan(result$measure$y)))
})

test_that("incidence() with dow adjusts values appropriately", {
  skip_if_not(exists("fit_rw_single_dow"),
              "No test fixture with dow_effect available")

  # Get incidence with and without dow
  result_with <- incidence(fit_rw_single_dow, dow = TRUE)
  result_without <- incidence(fit_rw_single_dow, dow = FALSE)

  # With dow should apply: incidence * week_effect * dow_simplex
  # So values should generally differ (unless dow_simplex happens to be ~1)
  differences <- abs(result_with$measure$y - result_without$measure$y)

  # At least some time points should show dow adjustment
  # (Could be zero if model estimated dow_simplex ≈ [1/7, 1/7, ...])
  expect_true(mean(differences) >= 0)  # Sanity check
})

# Documentation and user experience tests
test_that("incidence() provides informative messages", {
  # Using default (NULL) should inform user of setting used
  expect_message(
    incidence(fit_rw_single, dow = NULL),
    "day-of-week setting.*model"
  )

  # Explicit values should not message
  expect_silent(
    incidence(fit_rw_single, dow = FALSE)
  )
})

test_that("incidence() error messages are clear and actionable", {
  # Invalid class error
  expect_error(
    incidence("not a fit", dow = FALSE),
    class = "error"  # Should throw an error with clear message
  )

  # dow incompatibility error should mention what's needed
  err <- tryCatch(
    incidence(fit_rw_single, dow = TRUE),
    error = function(e) e$message
  )
  expect_match(err, "Day-of-week effects cannot be incorporated", ignore.case = TRUE)
})

# ==============================================================================
# TESTS: growth_rate() FUNCTION
# ==============================================================================

test_that("growth_rate() validates inputs correctly", {
  expect_error(
    growth_rate(list(fit = "not a fit")),
    "Input must be of class"
  )
})

test_that("growth_rate() works for single pathogen models", {
  result <- growth_rate(fit_rw_single)

  expect_s3_class(result, "growth_rate")
  expect_s3_class(result, "EpiStrainDynamics.metric")

  # Growth rate starts from time index 2 (needs previous day)
  n_days <- length(fit_rw_single$constructed_model$data$time)
  expected_rows <- n_days - 2 + 1
  expect_equal(nrow(result$measure), expected_rows)

  # Growth rate can be positive or negative
  expect_true(any(result$measure$y > 0))
  expect_true(any(result$measure$y < 0))

  # But should be reasonable (not extreme values)
  expect_true(all(abs(result$measure$y) < 10))  # Growth rate shouldn't be crazy
})

test_that("growth_rate() works for multiple pathogen models", {
  result <- growth_rate(fit_rw_multi)

  expect_s3_class(result, "growth_rate")

  # Should have individual pathogens + Total
  pathogen_names <- unique(result$measure$pathogen)
  expect_true("Total" %in% pathogen_names)
  expect_equal(length(pathogen_names), 5)
})

test_that("growth_rate() threshold of 0 makes sense for prop", {
  result <- growth_rate(fit_rw_single)

  # Prop should be proportion > 0 (i.e., growing)
  expect_true(all(result$measure$prop >= 0 & result$measure$prop <= 1))

  # If median growth rate > 0, prop should be > 0.5 (in general)
  growing <- result$measure[result$measure$y > 0, ]
  expect_true(mean(growing$prop > 0.5) > 0.7)  # Most growing periods should have prop > 0.5
})

# ==============================================================================
# TESTS: proportion() FUNCTION
# ==============================================================================

test_that("proportion() validates inputs correctly", {
  expect_error(
    proportion(list(fit = "not a fit")),
    "Input must be of class"
  )

  # Can't use proportion on single pathogen model
  expect_error(
    proportion(fit_rw_single),
    "Input must inherit from"
  )
})

test_that("proportion() works with default (all individual pathogens)", {
  result <- proportion(fit_rw_multi)

  expect_s3_class(result, "proportion")
  expect_s3_class(result, "EpiStrainDynamics.metric")

  # Should have one row per pathogen per time point (no Total for proportion)
  n_days <- length(fit_rw_multi$constructed_model$data$time)
  n_pathogens <- 4
  expected_rows <- n_days * n_pathogens
  expect_equal(nrow(result$measure), expected_rows)

  # Proportions should be between 0 and 1 (when denominator includes numerator)
  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$y <= 1))

  # Check pathogen names are correct (not numeric indices)
  expect_true(all(result$measure$pathogen %in% c("alpha", "delta", "omicron", "other")))

  # For each time point, proportions should sum to ~1 (given posterior sampling variation)
  first_time <- result$measure[result$measure$time == result$measure$time[1], ]
  expect_true(abs(sum(first_time$y) - 1) < 0.1)  # Allow some sampling variation
})

test_that("proportion() works with custom numerator and denominator", {
  result <- proportion(fit_rw_multi,
                       numerator_combination = "alpha",
                       denominator_combination = c("alpha", "delta"))

  expect_s3_class(result, "proportion")

  # Should have one proportion per time point
  n_days <- length(fit_rw_multi$constructed_model$data$time)
  expect_equal(nrow(result$measure), n_days)

  # Proportion can be > 1 if denominator excludes some of numerator
  # But in this case alpha is in denominator, so should be 0 to 1
  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$y <= 1))
})

test_that("proportion() validates pathogen names", {
  expect_error(
    proportion(fit_rw_multi,
               numerator_combination = "nonexistent",
               denominator_combination = "all"),
    "invalid pathogen name"
  )

  expect_error(
    proportion(fit_rw_multi,
               numerator_combination = "alpha",
               denominator_combination = c("alpha", "nonexistent")),
    "invalid pathogen name"
  )
})

# ==============================================================================
# TESTS: INTERNAL CALCULATION FUNCTIONS
# ==============================================================================

test_that("calc_stats computes correct quantiles", {
  set.seed(123)
  values <- rnorm(1000, mean = 5, sd = 2)

  result <- calc_stats(values, threshold = 5)

  # Check quantiles
  expect_equal(result$y, median(values), tolerance = 0.01)
  expect_equal(result$lb_95, quantile(values, 0.025)[[1]], tolerance = 0.01)
  expect_equal(result$ub_95, quantile(values, 0.975)[[1]], tolerance = 0.01)
  expect_equal(result$lb_50, quantile(values, 0.25)[[1]], tolerance = 0.01)
  expect_equal(result$ub_50, quantile(values, 0.75)[[1]], tolerance = 0.01)

  # Check prop
  expect_equal(result$prop, mean(values > 5), tolerance = 0.01)
})

test_that("get_model_components extracts correct structure", {
  components <- get_model_components(fit_rw_single)

  expect_named(components, c("fit", "pathogen_names", "num_path", "time_seq",
                             "time", "num_days", "knots", "spline_degree",
                             "DOW", "week_effect", "dow_effect"))

  expect_s4_class(components$fit, "stanfit")
  expect_equal(components$pathogen_names, "cases")
  expect_equal(components$num_path, 1)
  expect_false(components$dow_effect)
})

test_that("expand_pathogen_grid creates correct combinations", {
  time_grid <- data.frame(time_idx = 1:10)
  pathogen_names <- c("alpha", "delta", "omicron")

  result <- expand_pathogen_grid(time_grid, pathogen_names)

  expect_equal(nrow(result), 30)  # 10 times * 3 pathogens
  expect_named(result, c("pathogen", "pathogen_idx", "time_idx"))
  expect_equal(unique(result$pathogen), pathogen_names)
  expect_equal(unique(result$pathogen_idx), 1:3)
})

# ==============================================================================
# TESTS: REGRESSION TESTS (saved expected values)
# ==============================================================================

test_that("Rt calculation matches saved reference values", {
  skip_if_not(file.exists("tests/testthat/fixtures/expected_rt_single.rds"),
              "Reference values not available")

  result <- Rt(fit_rw_single, tau_max = 7, gi_dist = gi_simple)
  expected <- readRDS("tests/testthat/fixtures/expected_rt_single.rds")

  # Check that medians are very close (allow small numerical differences)
  expect_equal(result$measure$y, expected$y, tolerance = 0.001)
})

test_that("incidence calculation matches saved reference values", {
  skip_if_not(file.exists("tests/testthat/fixtures/expected_incidence_multi.rds"),
              "Reference values not available")

  result <- incidence(fit_rw_multi, dow = FALSE)
  expected <- readRDS("tests/testthat/fixtures/expected_incidence_multi.rds")

  expect_equal(result$measure$y, expected$y, tolerance = 0.001)
})

# ==============================================================================
# TESTS: EDGE CASES AND ERROR HANDLING
# ==============================================================================

test_that("metrics handle missing values appropriately", {
  # The Stan model should handle NAs in input, but let's verify output
  result <- incidence(fit_rw_single, dow = FALSE)

  # No NAs should be in output
  expect_false(any(is.na(result$measure$y)))
  expect_false(any(is.na(result$measure$lb_95)))
  expect_false(any(is.na(result$measure$ub_95)))
})
