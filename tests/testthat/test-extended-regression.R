# Regression tests for fitted models using pre-fitted fixtures downloaded
# via piggyback (handled by helper-fixtures.R).
#
#' @srrstats {G5.10} Extended tests included and run by default
#' @srrstats {G5.11, G5.11a} Extended fixtures downloaded by piggybank. Tests
#'   skip if download fails.
#' @srrstats {G5.12} Requirements described in tests/README.md

# ===========================================================================
# Rt estimation
# ===========================================================================

test_that("Rt returns expected structure from random walk single-pathogen model", {
  skip_if_no_extended_fixtures()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Rt point estimates are finite and non-negative", {
  skip_if_no_extended_fixtures()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  expect_true(all(is.finite(result$measure$y)))
  expect_true(all(result$measure$y >= 0))
})

test_that("Rt credible intervals are ordered correctly", {
  skip_if_no_extended_fixtures()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_rw_single, tau_max = 7, gi_dist = gi_simple)

  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$lb_50 <= result$measure$y))
  expect_true(all(result$measure$y <= result$measure$ub_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

# ===========================================================================
# Incidence
# ===========================================================================

test_that("Incidence returns expected structure from random walk multi-pathogen model", {
  skip_if_no_extended_fixtures()

  result <- incidence(extended_fixtures$fit_rw_multi, dow = FALSE)

  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Incidence estimates are finite and non-negative", {
  skip_if_no_extended_fixtures()

  result <- incidence(extended_fixtures$fit_rw_multi, dow = FALSE)

  expect_true(all(is.finite(result$measure$y)))
  expect_true(all(result$measure$y >= 0))
})

test_that("Incidence credible intervals are ordered correctly", {
  skip_if_no_extended_fixtures()

  result <- incidence(extended_fixtures$fit_rw_multi, dow = FALSE)

  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$lb_50 <= result$measure$y))
  expect_true(all(result$measure$y <= result$measure$ub_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

# ===========================================================================
# Growth rate
# ===========================================================================

test_that("Growth rate returns expected structure from random walk single-pathogen model", {
  skip_if_no_extended_fixtures()

  result <- growth_rate(extended_fixtures$fit_rw_single)

  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Growth rate estimates are finite", {
  skip_if_no_extended_fixtures()

  result <- growth_rate(extended_fixtures$fit_rw_single)

  expect_true(all(is.finite(result$measure$y)))
})

test_that("Growth rate credible intervals are ordered correctly", {
  skip_if_no_extended_fixtures()

  result <- growth_rate(extended_fixtures$fit_rw_single)

  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$lb_50 <= result$measure$y))
  expect_true(all(result$measure$y <= result$measure$ub_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

# ===========================================================================
# Proportion
# ===========================================================================

test_that("Proportion returns expected structure from random walk multi-pathogen model", {
  skip_if_no_extended_fixtures()

  result <- proportion(extended_fixtures$fit_rw_multi)

  expect_true(!is.null(result$measure))
  expect_true(nrow(result$measure) > 0)
})

test_that("Proportion estimates are bounded between 0 and 1", {
  skip_if_no_extended_fixtures()

  result <- proportion(extended_fixtures$fit_rw_multi)

  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$y <= 1))
})

test_that("Proportion credible intervals are ordered and bounded", {
  skip_if_no_extended_fixtures()

  result <- proportion(extended_fixtures$fit_rw_multi)

  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$lb_50 <= result$measure$y))
  expect_true(all(result$measure$y <= result$measure$ub_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
  expect_true(all(result$measure$lb_95 >= 0))
  expect_true(all(result$measure$ub_95 <= 1))
})

# ===========================================================================
# P-spline models produce comparable output structure to random walk
# ===========================================================================

test_that("P-spline single-pathogen model produces valid Rt estimates", {
  skip_if_no_extended_fixtures()

  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))

  result <- Rt(extended_fixtures$fit_ps_single, tau_max = 7, gi_dist = gi_simple)

  expect_true(!is.null(result$measure))
  expect_true(all(is.finite(result$measure$y)))
  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

test_that("P-spline multi-pathogen model produces valid proportion estimates", {
  skip_if_no_extended_fixtures()

  result <- proportion(extended_fixtures$fit_ps_multi)

  expect_true(!is.null(result$measure))
  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$y <= 1))
  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

# ===========================================================================
# Subtyped pathogen structure
# ===========================================================================

test_that("Random walk subtyped model produces valid incidence estimates", {
  skip_if_no_extended_fixtures()

  result <- incidence(extended_fixtures$fit_rw_subtyped, dow = FALSE)

  expect_true(!is.null(result$measure))
  expect_true(all(is.finite(result$measure$y)))
  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

test_that("P-spline subtyped model produces valid proportion estimates", {
  skip_if_no_extended_fixtures()

  result <- proportion(extended_fixtures$fit_ps_subtyped)

  expect_true(!is.null(result$measure))
  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$y <= 1))
  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})
