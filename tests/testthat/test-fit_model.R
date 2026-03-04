# Tests for fit_model() function
# Focus: fitting process, parameter validation, output structure
# Note: Detailed metric calculations are tested in test-metrics.R
# Note: Verbosity and error handling (BS2.12-BS2.15) tested in test-verbosity.R
#' @srrstats {G5.4, G5.5} Correctness tests with fixed test data + fixed seed

# Prevent browser from opening during tests (rstan diagnostics)
options(browser = function(...) {})

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

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2,
      case_timeseries = 'cases',
      time = 'date'
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
# TESTS: OUTPUT STRUCTURE (using cached fitted models from setup)
# ==============================================================================

test_that("fit_model() returns correct structure for single pathogen RW", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # Use the pre-fitted model from setup
  fit <- fit_rw_single

  # Check output structure
  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "rw_single")
  expect_named(fit, c("fit", "constructed_model", "mcmc_info"))

  # Check fit object is stanfit
  expect_s4_class(fit$fit, "stanfit")

  # Check constructed_model is present
  expect_s3_class(fit$constructed_model, "EpiStrainDynamics.model")
})

test_that("fit_model() returns correct structure for single pathogen PS", {
  skip_if_not(exists("fit_ps_single"), "Cached fitted models not available")

  fit <- fit_ps_single

  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "ps_single")
  expect_s4_class(fit$fit, "stanfit")
})

test_that("fit_model() returns correct structure for multiple pathogen RW", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  fit <- fit_rw_multi

  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "rw")
  expect_s4_class(fit$fit, "stanfit")
})

test_that("fit_model() returns correct structure for multiple pathogen PS", {
  skip_if_not(exists("fit_ps_multi"), "Cached fitted models not available")

  fit <- fit_ps_multi

  expect_s3_class(fit, "EpiStrainDynamics.fit")
  expect_s3_class(fit, "ps")
  expect_s4_class(fit$fit, "stanfit")
})

# ==============================================================================
# TESTS: MCMC PARAMETERS ARE RESPECTED
# ==============================================================================

test_that("fit_model() respects n_chain parameter", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # The cached model was fit with n_chain = 1
  expect_equal(fit_rw_single$fit@sim$chains, 1)
})

test_that("fit_model() respects n_iter parameter", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # The cached model was fit with n_iter = 500, n_warmup = 250 (default)
  expect_equal(fit_rw_single$fit@sim$iter, 500)
  expect_equal(fit_rw_single$fit@sim$warmup, 250)
})

# ==============================================================================
# TESTS: POSTERIOR SAMPLES HAVE CORRECT DIMENSIONS
# ==============================================================================

test_that("posterior samples have correct dimensions for single pathogen", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  post <- rstan::extract(fit_rw_single$fit)
  n_days <- nrow(fit_rw_single$constructed_model$validated_tsbl)

  # Check 'a' dimension: [samples, time]
  expect_equal(length(dim(post$a)), 2)
  expect_equal(ncol(post$a), n_days)

  # Check we have the right number of posterior samples
  # (n_iter - n_warmup) * n_chain = (500 - 250) * 1 = 250
  expect_equal(nrow(post$a), 250)
})

test_that("posterior samples have correct dimensions for multiple pathogens", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  post <- rstan::extract(fit_rw_multi$fit)
  n_days <- nrow(fit_rw_multi$constructed_model$validated_tsbl)

  # Check 'a' dimension: [samples, pathogens, time]
  expect_equal(length(dim(post$a)), 3)
  expect_equal(dim(post$a)[2], 4)        # n_pathogens
  expect_equal(dim(post$a)[3], n_days)   # n_days
  expect_equal(dim(post$a)[1], 250)      # posterior samples
})
