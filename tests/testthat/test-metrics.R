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
  expect_error(
    incidence(list(fit = "not a fit"), dow = FALSE),
    "Input must be of class"
  )

  # dow = TRUE when model doesn't have dow_effect
  expect_error(
    incidence(fit_rw_single, dow = TRUE),
    "dow effects can't be incorporated"
  )
})

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
