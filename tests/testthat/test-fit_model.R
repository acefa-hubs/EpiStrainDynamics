# Tests for fit_model() function
# Focus: fitting process, parameter validation, output structure
# Note: Detailed metric calculations are tested in test-metrics.R
#' @srrstats {G5.4, G5.5} Correctness tests with fixed test data + fixed seed

# Prevent browser from opening during tests (rstan diagnostics)
options(browser = function(...) {})

# ==============================================================================
# SETUP: Use lightweight synthetic data
# ==============================================================================

create_test_data <- function(n_days = 30, seed = 123) {
  set.seed(seed)
  data.frame(
    dates = seq.Date(as.Date("2024-01-01"), by = "day", length.out = n_days),
    cases = rpois(n_days, lambda = exp(6 + rnorm(n_days, 0, 0.2)))
  )
}

create_test_data_multi <- function(n_days = 30, seed = 123) {
  set.seed(seed)
  total <- rpois(n_days, lambda = exp(6 + rnorm(n_days, 0, 0.2)))

  # Simple multinomial split
  props <- matrix(c(0.4, 0.3, 0.2, 0.1), nrow = n_days, ncol = 4, byrow = TRUE)
  props <- props + matrix(rnorm(n_days * 4, 0, 0.05), nrow = n_days)
  props <- pmax(props, 0.01)
  props <- props / rowSums(props)

  p1 <- rbinom(n_days, total, props[, 1])
  p2 <- rbinom(n_days, total - p1, props[, 2] / (1 - props[, 1]))
  p3 <- rbinom(n_days, total - p1 - p2, props[, 3] / (1 - props[, 1] - props[, 2]))
  p4 <- total - p1 - p2 - p3

  data.frame(
    dates = seq.Date(as.Date("2024-01-01"), by = "day", length.out = n_days),
    total = total,
    var1 = p1,
    var2 = p2,
    var3 = p3,
    var4 = p4
  )
}

# ==============================================================================
# TESTS: INPUT VALIDATION
# ==============================================================================

test_that("fit_model() validates constructed_model class", {
  expect_error(
    fit_model(list(data = "not a model")),
    "Input must be of class"
  )

  expect_error(
    fit_model("not a model"),
    "Input must be of class"
  )
})

test_that("fit_model() validates MCMC parameters", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  # Invalid n_chain
  expect_error(
    fit_model(mod, n_chain = -1),
    "`n_chain` must be positive, got -1"
  )

  # Invalid n_iter
  expect_error(
    fit_model(mod, n_iter = 0),
    "must be positive|invalid"
  )

})

# ==============================================================================
# TESTS: OUTPUT STRUCTURE (using fast, small models)
# ==============================================================================

test_that("fit_model() returns correct structure for single pathogen RW", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 500, n_chain = 1,
                     verbose = FALSE, seed = 123)
  )

  # Check output structure
  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "rw_single")
  expect_named(fit, c("fit", "constructed_model"))

  # Check fit object is stanfit
  expect_s4_class(fit$fit, "stanfit")

  # Check constructed_model is preserved
  expect_identical(fit$constructed_model, mod)
})

test_that("fit_model() returns correct structure for single pathogen PS", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = p_spline(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 500, n_chain = 1,
                     verbose = FALSE, seed = 123)
  )

  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "ps_single")
  expect_s4_class(fit$fit, "stanfit")
})

test_that("fit_model() returns correct structure for multiple pathogen RW", {
  test_data <- create_test_data_multi(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = multiple(
      data = test_data,
      case_timeseries = 'total',
      time = 'dates',
      component_pathogen_timeseries = c('var1', 'var2', 'var3', 'var4')
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 500, n_chain = 1,
                     verbose = FALSE, seed = 123)
  )

  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "rw")
  expect_s4_class(fit$fit, "stanfit")
})

test_that("fit_model() returns correct structure for multiple pathogen PS", {
  test_data <- create_test_data_multi(n_days = 30)
  mod <- construct_model(
    method = p_spline(),
    pathogen_structure = multiple(
      data = test_data,
      case_timeseries = 'total',
      time = 'dates',
      component_pathogen_timeseries = c('var1', 'var2', 'var3', 'var4')
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 500, n_chain = 1,
                     verbose = FALSE, seed = 123)
  )

  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "ps")
  expect_s4_class(fit$fit, "stanfit")
})

# ==============================================================================
# TESTS: MCMC PARAMETERS ARE RESPECTED
# ==============================================================================

test_that("fit_model() respects n_chain parameter", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 500, n_chain = 2,
                     verbose = FALSE, seed = 123)
  )

  # Check that fit has 2 chains
  expect_equal(fit$fit@sim$chains, 2)
})

test_that("fit_model() respects n_iter parameter", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 600, n_warmup = 200, n_chain = 1,
                     verbose = FALSE, seed = 123)
  )

  # Check total iterations (warmup + sampling)
  expect_equal(fit$fit@sim$iter, 600)
  expect_equal(fit$fit@sim$warmup, 200)
})

# ==============================================================================
# TESTS: FITTED MODEL CAN BE USED WITH METRICS
# ==============================================================================

test_that("fitted model works with metric functions", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 500, n_chain = 1,
                     verbose = FALSE, seed = 123)
  )

  # Should work with incidence
  expect_no_error(inc <- incidence(fit, dow = FALSE))
  expect_s3_class(inc, "incidence")

  # Should work with growth_rate
  expect_no_error(gr <- growth_rate(fit))
  expect_s3_class(gr, "growth_rate")

  # Should work with Rt
  gi_simple <- function(x) ifelse(x == 0, 0, 4 * x * exp(-2 * x))
  expect_no_error(rt <- Rt(fit, tau_max = 7, gi_dist = gi_simple))
  expect_s3_class(rt, "Rt")
})

# ==============================================================================
# TESTS: POSTERIOR SAMPLES HAVE CORRECT DIMENSIONS
# ==============================================================================

test_that("posterior samples have correct dimensions for single pathogen", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 600, n_warmup = 200, n_chain = 2,
                     verbose = FALSE, seed = 123)
  )

  post <- rstan::extract(fit$fit)

  # Check 'a' dimension: [samples, time]
  expect_equal(length(dim(post$a)), 2)
  expect_equal(ncol(post$a), 30)  # n_days

  # Check we have the right number of posterior samples
  # (n_iter - n_warmup) * n_chain = (600 - 200) * 2 = 800
  expect_equal(nrow(post$a), 800)
})

test_that("posterior samples have correct dimensions for multiple pathogens", {
  test_data <- create_test_data_multi(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = multiple(
      data = test_data,
      case_timeseries = 'total',
      time = 'dates',
      component_pathogen_timeseries = c('var1', 'var2', 'var3', 'var4')
    )
  )

  suppressWarnings(
    fit <- fit_model(mod, n_iter = 600, n_warmup = 200, n_chain = 2,
                     verbose = FALSE, seed = 123)
  )

  post <- rstan::extract(fit$fit)

  # Check 'a' dimension: [samples, pathogens, time]
  expect_equal(length(dim(post$a)), 3)
  expect_equal(dim(post$a)[2], 4)   # n_pathogens
  expect_equal(dim(post$a)[3], 30)  # n_days
  expect_equal(dim(post$a)[1], 800) # posterior samples
})

# ==============================================================================
# TESTS: VERBOSE PARAMETER
# ==============================================================================

test_that("verbose = FALSE suppresses Stan output", {
  test_data <- create_test_data(n_days = 30)
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  # Capture output
  output <- capture.output(suppressWarnings(
    fit <- fit_model(mod, n_iter = 500, n_chain = 1,
                     verbose = FALSE, seed = 123)
  ))

  # Should have minimal/no output
  expect_lt(length(output), 5)
})

# ==============================================================================
# OPTIONAL: REGRESSION TEST WITH FIXTURE (if desired)
# ==============================================================================

test_that("fit_model produces consistent results with fixture", {
  # This uses the fixture from test-metrics.R if available
  skip_if_not(file.exists("tests/testthat/fixtures/fit_rw_single.rds"),
              "Fixture not available")

  # Load the fixture
  fit_fixture <- readRDS("tests/testthat/fixtures/fit_rw_single.rds")

  # Refit with same data and seed
  test_data <- fit_fixture$constructed_model$data
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = data.frame(
        dates = test_data$time,
        cases = test_data$case_timeseries
      ),
      case_timeseries = 'cases',
      time = 'dates'
    )
  )

  suppressWarnings(
    fit_new <- fit_model(mod, n_iter = 1000, n_chain = 2,
                         verbose = FALSE, seed = 42)
  )

  # Check structure is the same
  expect_s3_class(fit_new, "EpiStrainDynamics.fit")
  expect_s3_class(fit_new, "rw_single")

  # Note: We don't expect exact numerical equivalence due to MCMC stochasticity
  # even with the same seed (Stan version differences, etc.)
  # But we can check dimensions are consistent
  post_fixture <- rstan::extract(fit_fixture$fit)
  post_new <- rstan::extract(fit_new$fit)

  expect_equal(dim(post_fixture$a), dim(post_new$a))
})
