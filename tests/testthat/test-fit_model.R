test_that("fit_model parameters are passed correctly to Stan", {
  skip_if_not_installed("rstan")

  # Capture parameters passed to rstan::sampling
  captured_params <- list()

  # Create custom mock in test environment
  test_env$mock_sampling_capture <- function(object, data, iter, warmup, chains) {
    captured_params$iter <<- iter
    captured_params$warmup <<- warmup
    captured_params$chains <<- chains
    return(test_env$mock_stan_fit)
  }

  # Set up mocks with parameter capture
  setup_stan_mocks()
  assignInNamespace("sampling", test_env$mock_sampling_capture, ns = "rstan")
  assign("stanmodels", test_env$mock_stanmodels, envir = .GlobalEnv)
  on.exit(teardown_mocks())

  models <- create_test_models()
  fit_model(models$rw_single, iter = 500, warmup = 200, chains = 4)

  expect_equal(captured_params$iter, 500)
  expect_equal(captured_params$warmup, 200)
  expect_equal(captured_params$chains, 4)
})

test_that("standata structure is correct for different model types", {
  skip_if_not_installed("rstan")

  captured_data <- list()

  # Create custom mock in test environment
  test_env$mock_sampling_data <- function(object, data, iter, warmup, chains) {
    captured_data <<- data
    return(test_env$mock_stan_fit)
  }

  setup_stan_mocks()
  assignInNamespace("sampling", test_env$mock_sampling_data, ns = "rstan")
  assign("stanmodels", test_env$mock_stanmodels, envir = .GlobalEnv)
  on.exit(teardown_mocks())

  models <- create_test_models()
  expected_lengths <- get_expected_data_lengths()

  # Test single model standata
  fit_model(models$rw_single, iter = 10, warmup = 5, chains = 1)

  expect_true("num_data" %in% names(captured_data))
  expect_true("Y" %in% names(captured_data))
  expect_true("week_effect" %in% names(captured_data))
  expect_true("DOW" %in% names(captured_data))
  expect_equal(length(captured_data$Y), expected_lengths$sarscov2_length)

  # Test multiple model standata
  fit_model(models$rw_multiple, iter = 10, warmup = 5, chains = 1)

  expect_true("num_path" %in% names(captured_data))
  expect_true("P" %in% names(captured_data))
  expect_true("cov_structure" %in% names(captured_data))
  expect_true("noise_structure" %in% names(captured_data))
  expect_equal(captured_data$num_path, expected_lengths$sarscov2_pathogen_count)

  # Test subtyped model standata
  fit_model(models$rw_subtyped, iter = 10, warmup = 5, chains = 1)

  expect_true("P1" %in% names(captured_data))
  expect_true("P2" %in% names(captured_data))
  expect_equal(length(captured_data$Y), expected_lengths$influenza_length)
})

test_that("model-specific standata validation", {
  skip_if_not_installed("rstan")

  all_captured_data <- list()

  # Create custom mock in test environment
  test_env$mock_sampling_data_collector <- function(object, data, iter, warmup, chains) {
    # Store data with a unique key
    key <- paste0("call_", length(all_captured_data) + 1)
    all_captured_data[[key]] <<- data
    return(test_env$mock_stan_fit)
  }

  setup_stan_mocks()
  assignInNamespace("sampling", test_env$mock_sampling_data_collector, ns = "rstan")
  assign("stanmodels", test_env$mock_stanmodels, envir = .GlobalEnv)
  on.exit(teardown_mocks())

  models <- create_test_models()
  expected_lengths <- get_expected_data_lengths()
  expected_names <- get_expected_pathogen_names()

  # Test single models have correct structure
  fit_model(models$rw_single, iter = 10, warmup = 5, chains = 1)
  single_data <- all_captured_data[[length(all_captured_data)]]
  expect_false("num_path" %in% names(single_data))
  expect_false("P" %in% names(single_data))
  expect_equal(length(single_data$Y), expected_lengths$sarscov2_length)

  # Test P-spline models have spline parameters
  fit_model(models$ps_single, iter = 10, warmup = 5, chains = 1)
  ps_data <- all_captured_data[[length(all_captured_data)]]
  expect_true("X" %in% names(ps_data))
  expect_true("num_knots" %in% names(ps_data))
  expect_equal(length(ps_data$Y), expected_lengths$sarscov2_length)

  # Test multiple models have pathogen data
  fit_model(models$rw_multiple, iter = 10, warmup = 5, chains = 1)
  multiple_data <- all_captured_data[[length(all_captured_data)]]
  expect_true("num_path" %in% names(multiple_data))
  expect_true("P" %in% names(multiple_data))
  expect_equal(multiple_data$num_path, expected_lengths$sarscov2_pathogen_count)
  expect_equal(length(multiple_data$Y), expected_lengths$sarscov2_length)

  # Test subtyped models have P1 and P2
  fit_model(models$rw_subtyped, iter = 10, warmup = 5, chains = 1)
  subtyped_data <- all_captured_data[[length(all_captured_data)]]
  expect_true("P1" %in% names(subtyped_data))
  expect_true("P2" %in% names(subtyped_data))
  expect_false("P" %in% names(subtyped_data)) # Should not have single P matrix
  expect_equal(length(subtyped_data$Y), expected_lengths$influenza_length)
})
