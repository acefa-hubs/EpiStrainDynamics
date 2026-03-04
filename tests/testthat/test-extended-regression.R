# Regression tests for fitted models. These tests are slow because they
# depend on model fitting (handled by helper-extended-fixtures.R).
#
# To run:
#   EPISTRAINDYNAMICS_EXTENDED_TESTS=true Rscript -e 'devtools::test()'
#
#' @srrstats {G5.10} Extended tests included and switched on by flag
#' @srrstats {G5.11, G5.11a} All data is package data (sarscov2, influenza),
#'   fully reproducible without downloading external assets.
#' @srrstats {G5.12} Requirements described in tests/README.md
#'
# ===========================================================================
# Rt estimation
# ===========================================================================

test_that("Rt returns expected structure from random walk single-pathogen model", {
  skip_extended()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  # Structure checks
  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Rt point estimates are finite and non-negative", {
  skip_extended()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  # All point estimates should be finite
  expect_true(all(is.finite(result$measure$median)))

  # Rt should be non-negative
  expect_true(all(result$measure$median >= 0))
})

test_that("Rt credible intervals are ordered correctly", {
  skip_extended()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  # Lower bound <= median <= upper bound
  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
})

# ===========================================================================
# Incidence
# ===========================================================================

test_that("Incidence returns expected structure from random walk multi-pathogen model", {
  skip_extended()

  result <- incidence(extended_fixtures$fit_rw_multi, dow = FALSE)

  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Incidence estimates are finite and non-negative", {
  skip_extended()

  result <- incidence(extended_fixtures$fit_rw_multi, dow = FALSE)

  expect_true(all(is.finite(result$measure$median)))
  expect_true(all(result$measure$median >= 0))
})

test_that("Incidence credible intervals are ordered correctly", {
  skip_extended()

  result <- incidence(extended_fixtures$fit_rw_multi, dow = FALSE)

  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
})

# ===========================================================================
# Growth rate
# ===========================================================================

test_that("Growth rate returns expected structure from random walk single-pathogen model", {
  skip_extended()

  result <- growth_rate(extended_fixtures$fit_rw_single)

  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Growth rate estimates are finite", {
  skip_extended()

  result <- growth_rate(extended_fixtures$fit_rw_single)

  # Growth rates can be negative, but should be finite
  expect_true(all(is.finite(result$measure$median)))
})

test_that("Growth rate credible intervals are ordered correctly", {
  skip_extended()

  result <- growth_rate(extended_fixtures$fit_rw_single)

  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
})

# ===========================================================================
# Proportion
# ===========================================================================

test_that("Proportion returns expected structure from random walk multi-pathogen model", {
  skip_extended()

  result <- proportion(extended_fixtures$fit_rw_multi)

  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Proportion estimates are bounded between 0 and 1", {
  skip_extended()

  result <- proportion(extended_fixtures$fit_rw_multi)

  expect_true(all(result$measure$median >= 0))
  expect_true(all(result$measure$median <= 1))
})

test_that("Proportion credible intervals are ordered and bounded", {
  skip_extended()

  result <- proportion(extended_fixtures$fit_rw_multi)

  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
  expect_true(all(result$measure$lower >= 0))
  expect_true(all(result$measure$upper <= 1))
})
#
# # ===========================================================================
# # Day-of-week model
# # ===========================================================================
#
# test_that("Day-of-week model produces finite, non-negative incidence estimates", {
#   skip_extended()
#
#   result <- incidence(extended_fixtures$fit_rw_single_dow, dow = TRUE)
#
#   expect_true(!is.null(result$measure))
#   expect_true(all(is.finite(result$measure$median)))
#   expect_true(all(result$measure$median >= 0))
# })

# ===========================================================================
# P-spline models produce comparable output structure to random walk
# ===========================================================================

test_that("P-spline single-pathogen model produces valid Rt estimates", {
  skip_extended()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_ps_single, tau_max = 7, gi_dist = gi_simple)

  expect_true(!is.null(result$measure))
  expect_true(all(is.finite(result$measure$median)))
  expect_true(all(result$measure$median >= 0))
  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
})

test_that("P-spline multi-pathogen model produces valid proportion estimates", {
  skip_extended()

  result <- proportion(extended_fixtures$fit_ps_multi)

  expect_true(!is.null(result$measure))
  expect_true(all(result$measure$median >= 0))
  expect_true(all(result$measure$median <= 1))
  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
})

# ===========================================================================
# Subtyped pathogen structure
# ===========================================================================

test_that("Random walk subtyped model produces valid incidence estimates", {
  skip_extended()

  result <- incidence(extended_fixtures$fit_rw_subtyped, dow = FALSE)

  expect_true(!is.null(result$measure))
  expect_true(all(is.finite(result$measure$median)))
  expect_true(all(result$measure$median >= 0))
  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
})

test_that("P-spline subtyped model produces valid proportion estimates", {
  skip_extended()

  result <- proportion(extended_fixtures$fit_ps_subtyped)

  expect_true(!is.null(result$measure))
  expect_true(all(result$measure$median >= 0))
  expect_true(all(result$measure$median <= 1))
  expect_true(all(result$measure$lower <= result$measure$median))
  expect_true(all(result$measure$median <= result$measure$upper))
})
