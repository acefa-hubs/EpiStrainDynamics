#' @srrstats {G5.2, G5.2b} error and warning behaviour explicitly tested,
#'   including conditions that trigger these messages and comparing with
#'   expected results
#' @srrstats {G3.0} comparisons made between appropriate values

# ==============================================================================
# TESTS: validate_n_chain()
# ==============================================================================

test_that("validate_n_chain() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_n_chain(1))
  expect_silent(EpiStrainDynamics:::validate_n_chain(4))
  expect_silent(EpiStrainDynamics:::validate_n_chain(1L))
})

test_that("validate_n_chain() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_n_chain("4"), "must be numeric")
  expect_error(EpiStrainDynamics:::validate_n_chain(TRUE), "must be numeric")
})

test_that("validate_n_chain() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_n_chain(c(1, 2)), "must be a single value")
})

test_that("validate_n_chain() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_n_chain(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_n_chain(NA_real_), "must be a finite number")
})

test_that("validate_n_chain() rejects non-whole numbers", {
  expect_error(EpiStrainDynamics:::validate_n_chain(1.5), "must be a whole number")
})

test_that("validate_n_chain() rejects non-positive values", {
  expect_error(EpiStrainDynamics:::validate_n_chain(0), "must be positive")
  expect_error(EpiStrainDynamics:::validate_n_chain(-1), "must be positive")
})

# ==============================================================================
# TESTS: validate_n_iter()
# ==============================================================================

test_that("validate_n_iter() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_n_iter(100))
  expect_silent(EpiStrainDynamics:::validate_n_iter(2000L))
})

test_that("validate_n_iter() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_n_iter("100"), "must be numeric")
})

test_that("validate_n_iter() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_n_iter(c(100, 200)), "must be a single value")
})

test_that("validate_n_iter() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_n_iter(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_n_iter(NA_real_), "must be a finite number")
})

test_that("validate_n_iter() rejects non-whole numbers", {
  expect_error(EpiStrainDynamics:::validate_n_iter(100.5), "must be a whole number")
})

test_that("validate_n_iter() rejects non-positive values", {
  expect_error(EpiStrainDynamics:::validate_n_iter(0), "must be positive")
  expect_error(EpiStrainDynamics:::validate_n_iter(-100), "must be positive")
})

# ==============================================================================
# TESTS: validate_n_warmup()
# ==============================================================================

test_that("validate_n_warmup() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_n_warmup(0))
  expect_silent(EpiStrainDynamics:::validate_n_warmup(500))
  expect_silent(EpiStrainDynamics:::validate_n_warmup(500L))
})

test_that("validate_n_warmup() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_n_warmup("500"), "must be numeric")
})

test_that("validate_n_warmup() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_n_warmup(c(100, 200)), "must be a single value")
})

test_that("validate_n_warmup() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_n_warmup(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_n_warmup(NA_real_), "must be a finite number")
})

test_that("validate_n_warmup() rejects non-whole numbers", {
  expect_error(EpiStrainDynamics:::validate_n_warmup(100.5), "must be a whole number")
})

test_that("validate_n_warmup() rejects negative values", {
  expect_error(EpiStrainDynamics:::validate_n_warmup(-1), "cannot be negative")
})

# ==============================================================================
# TESTS: validate_thin()
# ==============================================================================

test_that("validate_thin() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_thin(1))
  expect_silent(EpiStrainDynamics:::validate_thin(5L))
})

test_that("validate_thin() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_thin("1"), "must be numeric")
})

test_that("validate_thin() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_thin(c(1, 2)), "must be a single value")
})

test_that("validate_thin() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_thin(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_thin(NA_real_), "must be a finite number")
})

test_that("validate_thin() rejects non-whole numbers", {
  expect_error(EpiStrainDynamics:::validate_thin(1.5), "must be a whole number")
})

test_that("validate_thin() rejects non-positive values", {
  expect_error(EpiStrainDynamics:::validate_thin(0), "must be positive")
  expect_error(EpiStrainDynamics:::validate_thin(-1), "must be positive")
})

# ==============================================================================
# TESTS: validate_adapt_delta()
# ==============================================================================

test_that("validate_adapt_delta() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_adapt_delta(0.9))
  expect_silent(EpiStrainDynamics:::validate_adapt_delta(0.5))
  expect_silent(EpiStrainDynamics:::validate_adapt_delta(0.01))
  expect_silent(EpiStrainDynamics:::validate_adapt_delta(0.99))
})

test_that("validate_adapt_delta() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_adapt_delta("0.9"), "must be numeric")
})

test_that("validate_adapt_delta() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_adapt_delta(c(0.9, 0.8)), "must be a single value")
})

test_that("validate_adapt_delta() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_adapt_delta(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_adapt_delta(NA_real_), "must be a finite number")
})

test_that("validate_adapt_delta() rejects out-of-range values", {
  expect_error(EpiStrainDynamics:::validate_adapt_delta(0), "must be between 0 and 1")
  expect_error(EpiStrainDynamics:::validate_adapt_delta(1), "must be between 0 and 1")
  expect_error(EpiStrainDynamics:::validate_adapt_delta(-0.1), "must be between 0 and 1")
  expect_error(EpiStrainDynamics:::validate_adapt_delta(1.1), "must be between 0 and 1")
})

# ==============================================================================
# TESTS: validate_seed()
# ==============================================================================

test_that("validate_seed() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_seed(NULL))
  expect_silent(EpiStrainDynamics:::validate_seed(123))
  expect_silent(EpiStrainDynamics:::validate_seed(1L))
})

test_that("validate_seed() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_seed("123"), "must be numeric or NULL")
})

test_that("validate_seed() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_seed(c(1, 2)), "must be a single value")
})

test_that("validate_seed() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_seed(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_seed(NA_real_), "must be a finite number")
})

test_that("validate_seed() rejects non-whole numbers", {
  expect_error(EpiStrainDynamics:::validate_seed(123.5), "must be a whole number")
})

test_that("validate_seed() rejects non-positive values", {
  expect_error(EpiStrainDynamics:::validate_seed(0), "must be positive")
  expect_error(EpiStrainDynamics:::validate_seed(-1), "must be positive")
})

# ==============================================================================
# TESTS: validate_verbose()
# ==============================================================================

test_that("validate_verbose() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_verbose(TRUE))
  expect_silent(EpiStrainDynamics:::validate_verbose(FALSE))
})

test_that("validate_verbose() rejects non-logical", {
  expect_error(EpiStrainDynamics:::validate_verbose(1), "must be logical")
  expect_error(EpiStrainDynamics:::validate_verbose("TRUE"), "must be logical")
})

test_that("validate_verbose() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_verbose(c(TRUE, FALSE)), "must be a single value")
})

test_that("validate_verbose() rejects NA", {
  expect_error(EpiStrainDynamics:::validate_verbose(NA), "cannot be NA")
})

# ==============================================================================
# TESTS: validate_multi_cores()
# ==============================================================================

test_that("validate_multi_cores() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_multi_cores(TRUE))
  expect_silent(EpiStrainDynamics:::validate_multi_cores(FALSE))
})

test_that("validate_multi_cores() rejects non-logical", {
  expect_error(EpiStrainDynamics:::validate_multi_cores(1), "must be logical")
  expect_error(EpiStrainDynamics:::validate_multi_cores("TRUE"), "must be logical")
})

test_that("validate_multi_cores() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_multi_cores(c(TRUE, FALSE)), "must be a single value")
})

test_that("validate_multi_cores() rejects NA", {
  expect_error(EpiStrainDynamics:::validate_multi_cores(NA), "cannot be NA")
})

# ==============================================================================
# TESTS: validate_suppress_warnings()
# ==============================================================================

test_that("validate_suppress_warnings() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_suppress_warnings(TRUE))
  expect_silent(EpiStrainDynamics:::validate_suppress_warnings(FALSE))
})

test_that("validate_suppress_warnings() rejects non-logical", {
  expect_error(EpiStrainDynamics:::validate_suppress_warnings(1), "must be logical")
  expect_error(EpiStrainDynamics:::validate_suppress_warnings("TRUE"), "must be logical")
})

test_that("validate_suppress_warnings() rejects multiple values", {
  expect_error(
    EpiStrainDynamics:::validate_suppress_warnings(c(TRUE, FALSE)),
    "must be a single value"
  )
})

test_that("validate_suppress_warnings() rejects NA", {
  expect_error(EpiStrainDynamics:::validate_suppress_warnings(NA), "cannot be NA")
})

# ==============================================================================
# TESTS: validate_rhat_threshold()
# ==============================================================================

test_that("validate_rhat_threshold() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_rhat_threshold(1.0))
  expect_silent(EpiStrainDynamics:::validate_rhat_threshold(1.01))
  expect_silent(EpiStrainDynamics:::validate_rhat_threshold(1.1))
})

test_that("validate_rhat_threshold() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_rhat_threshold("1.1"), "must be numeric")
})

test_that("validate_rhat_threshold() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_rhat_threshold(c(1.01, 1.1)), "must be a single value")
})

test_that("validate_rhat_threshold() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_rhat_threshold(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_rhat_threshold(NA_real_), "must be a finite number")
})

test_that("validate_rhat_threshold() rejects values below 1.0", {
  expect_error(
    EpiStrainDynamics:::validate_rhat_threshold(0.99),
    "must be greater than or equal to 1.0"
  )
  expect_error(
    EpiStrainDynamics:::validate_rhat_threshold(0.5),
    "must be greater than or equal to 1.0"
  )
})

test_that("validate_rhat_threshold() warns on unusually high values", {
  expect_warning(EpiStrainDynamics:::validate_rhat_threshold(5.0), "unusually high")
})

# ==============================================================================
# TESTS: validate_eff_sample_threshold()
# ==============================================================================

test_that("validate_eff_sample_threshold() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_eff_sample_threshold(100))
  expect_silent(EpiStrainDynamics:::validate_eff_sample_threshold(400))
  expect_silent(EpiStrainDynamics:::validate_eff_sample_threshold(1000))
})

test_that("validate_eff_sample_threshold() rejects non-numeric", {
  expect_error(EpiStrainDynamics:::validate_eff_sample_threshold("100"), "must be numeric")
})

test_that("validate_eff_sample_threshold() rejects multiple values", {
  expect_error(
    EpiStrainDynamics:::validate_eff_sample_threshold(c(100, 200)),
    "must be a single value"
  )
})

test_that("validate_eff_sample_threshold() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_eff_sample_threshold(Inf), "must be a finite number")
  expect_error(EpiStrainDynamics:::validate_eff_sample_threshold(NA_real_), "must be a finite number")
})

test_that("validate_eff_sample_threshold() rejects non-positive values", {
  expect_error(EpiStrainDynamics:::validate_eff_sample_threshold(0), "must be positive")
  expect_error(EpiStrainDynamics:::validate_eff_sample_threshold(-100), "must be positive")
})

test_that("validate_eff_sample_threshold() warns on very low values", {
  expect_warning(EpiStrainDynamics:::validate_eff_sample_threshold(10), "very low")
})

test_that("validate_eff_sample_threshold() warns on very high values", {
  expect_warning(EpiStrainDynamics:::validate_eff_sample_threshold(10000), "very high")
})

# ==============================================================================
# TESTS: validate_mcmc_params_collective()
# ==============================================================================

test_that("validate_mcmc_params_collective() checks n_warmup < n_iter", {
  expect_error(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 100, n_warmup = 100, n_chain = 1
    ),
    "n_warmup.*must be less than.*n_iter"
  )
  expect_error(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 100, n_warmup = 150, n_chain = 1
    ),
    "n_warmup.*must be less than.*n_iter"
  )
})

test_that("validate_mcmc_params_collective() warns about large n_iter", {
  expect_warning(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 25000, n_warmup = 1000, n_chain = 1
    ),
    "Large.*n_iter.*long computation"
  )
  expect_silent(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 2000, n_warmup = 1000, n_chain = 1
    )
  )
})

test_that("validate_mcmc_params_collective() warns when n_chain exceeds cores", {
  n_cores <- parallel::detectCores(logical = FALSE)
  if (!is.na(n_cores)) {
    expect_warning(
      EpiStrainDynamics:::validate_mcmc_params_collective(
        n_iter = 1000, n_warmup = 500, n_chain = n_cores + 2
      ),
      "exceeds available CPU cores"
    )
  }
})

test_that("validate_mcmc_params_collective() warns about excessive thinning", {
  expect_warning(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 100, n_warmup = 50, n_chain = 1, thin = 5
    ),
    "Very few effective samples.*thinning"
  )
})

test_that("validate_mcmc_params_collective() validates seed parameter", {
  expect_silent(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 1000, n_warmup = 500, n_chain = 1, seed = NULL
    )
  )
  expect_silent(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 1000, n_warmup = 500, n_chain = 1, seed = 123
    )
  )
  expect_error(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 1000, n_warmup = 500, n_chain = 1, seed = "abc"
    ),
    "seed.*must be a single numeric value"
  )
  expect_error(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 1000, n_warmup = 500, n_chain = 1, seed = c(1, 2)
    ),
    "seed.*must be a single numeric value"
  )
  expect_warning(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 1000, n_warmup = 500, n_chain = 1, seed = 123.5
    ),
    "seed.*not an integer"
  )
})

test_that("validate_mcmc_params_collective() respects suppress_warnings", {
  expect_silent(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 100, n_warmup = 50, n_chain = 1, thin = 5,
      suppress_warnings = TRUE
    )
  )
  expect_silent(
    EpiStrainDynamics:::validate_mcmc_params_collective(
      n_iter = 25000, n_warmup = 1000, n_chain = 1,
      suppress_warnings = TRUE
    )
  )
})

test_that("fit_model uses collective validation", {
  test_data <- data.frame(
    dates = seq.Date(as.Date("2024-01-01"), by = "day", length.out = 30),
    cases = rpois(30, lambda = 100)
  )
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = test_data,
      case_timeseries = "cases",
      time = "dates"
    )
  )
  expect_error(
    fit_model(mod, n_iter = 100, n_warmup = 100, n_chain = 1, verbose = FALSE),
    "n_warmup.*must be less than.*n_iter"
  )
})
