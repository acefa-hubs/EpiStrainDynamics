# Helper functions for creating test data and models using package datasets

# Create environment for storing mock objects
test_env <- new.env(parent = emptyenv())

# Helper function to check if package data exists
check_package_data <- function() {
  if (!exists("sarscov2") || !exists("influenza")) {
    skip("Package data (sarscov2, influenza) not available")
  }
}

# Helper function for validation testing using error expectations
validate_length_mismatch <- function(constructor_func, base_args, mismatch_args, expected_pattern = "length|same|match") {
  args <- modifyList(base_args, mismatch_args)
  expect_error(do.call(constructor_func, args), expected_pattern, ignore.case = TRUE)
}

# Helper function for testing method constructor structure
test_method_structure <- function(result, expected_method, expected_params = NULL) {
  # Test basic structure
  expect_type(result, "list")
  expect_s3_class(result, "EpiStrainDynamics.method")
  expect_equal(result$method, expected_method)
  expect_true("method" %in% names(result))
  expect_type(result$method, "character")

  # Test parameters if expected
  if (!is.null(expected_params)) {
    expect_true("model_params" %in% names(result))
    expect_type(result$model_params, "list")
    expect_named(result$model_params, names(expected_params))

    for (param_name in names(expected_params)) {
      expect_equal(result$model_params[[param_name]], expected_params[[param_name]])
    }
  } else {
    expect_length(result, 1)
  }
}

# Create standard test models using package datasets
create_test_models <- function() {
  check_package_data()

  list(
    rw_single = construct_model(
      method = random_walk(),
      pathogen_structure = single(
        case_timeseries = sarscov2$cases,
        time = sarscov2$date
      )
    ),

    ps_single = construct_model(
      method = p_spline(),
      pathogen_structure = single(
        case_timeseries = sarscov2$cases,
        time = sarscov2$date
      )
    ),

    rw_multiple = construct_model(
      method = random_walk(),
      pathogen_structure = multiple(
        case_timeseries = sarscov2$cases,
        time = sarscov2$date,
        component_pathogen_timeseries = list(
          alpha = sarscov2$alpha,
          delta = sarscov2$delta,
          omicron = sarscov2$omicron,
          other = sarscov2$other
        )
      )
    ),

    ps_multiple = construct_model(
      method = p_spline(),
      pathogen_structure = multiple(
        case_timeseries = sarscov2$cases,
        time = sarscov2$date,
        component_pathogen_timeseries = list(
          alpha = sarscov2$alpha,
          delta = sarscov2$delta,
          omicron = sarscov2$omicron,
          other = sarscov2$other
        )
      )
    ),

    rw_subtyped = construct_model(
      method = random_walk(),
      pathogen_structure = subtyped(
        case_timeseries = influenza$ili,
        time = influenza$week,
        influenzaA_unsubtyped_timeseries = influenza$inf_A,
        influenzaA_subtyped_timeseries = list(
          H3N2 = influenza$inf_H3N2,
          H1N1 = influenza$inf_H1N1
        ),
        other_pathogen_timeseries = list(
          influenzaB = influenza$inf_B,
          other = influenza$num_spec - influenza$inf_all
        )
      )
    ),

    ps_subtyped = construct_model(
      method = p_spline(),
      pathogen_structure = subtyped(
        case_timeseries = influenza$ili,
        time = influenza$week,
        influenzaA_unsubtyped_timeseries = influenza$inf_A,
        influenzaA_subtyped_timeseries = list(
          H3N2 = influenza$inf_H3N2,
          H1N1 = influenza$inf_H1N1
        ),
        other_pathogen_timeseries = list(
          influenzaB = influenza$inf_B,
          other = influenza$num_spec - influenza$inf_all
        )
      )
    )
  )
}

# Mock setup functions for Stan
setup_stan_mocks <- function() {
  # Create mock objects in the test environment
  test_env$mock_stanmodels <- list(
    rw_single = "mock_model",
    ps_single = "mock_model",
    rw_multiple = "mock_model",
    ps_multiple = "mock_model",
    rw_subtyped = "mock_model",
    ps_subtyped = "mock_model"
  )

  test_env$mock_stan_fit <- structure(list(), class = "stanfit")

  test_env$mock_rstan_sampling <- function(object, data, iter, warmup, chains) {
    return(test_env$mock_stan_fit)
  }

  # Store originals in test environment
  test_env$original_sampling <- NULL
  if (exists("sampling", envir = getNamespace("rstan"))) {
    test_env$original_sampling <- get("sampling", envir = getNamespace("rstan"))
  }

  test_env$original_stanmodels <- NULL
  if (exists("stanmodels", envir = .GlobalEnv)) {
    test_env$original_stanmodels <- get("stanmodels", envir = .GlobalEnv)
  }
}

setup_mocks <- function() {
  setup_stan_mocks()

  # Use the mock functions from test environment
  assignInNamespace("sampling", test_env$mock_rstan_sampling, ns = "rstan")
  assign("stanmodels", test_env$mock_stanmodels, envir = .GlobalEnv)
}

teardown_mocks <- function() {
  if (!is.null(test_env$original_sampling)) {
    assignInNamespace("sampling", test_env$original_sampling, ns = "rstan")
  }
  if (!is.null(test_env$original_stanmodels)) {
    assign("stanmodels", test_env$original_stanmodels, envir = .GlobalEnv)
  } else {
    if (exists("stanmodels", envir = .GlobalEnv)) {
      rm("stanmodels", envir = .GlobalEnv)
    }
  }
}

# Helper function to get expected data lengths from package datasets
get_expected_data_lengths <- function() {
  check_package_data()

  list(
    sarscov2_length = length(sarscov2$cases),
    influenza_length = length(influenza$ili),
    sarscov2_pathogen_count = 4,  # alpha, delta, omicron, other
    influenza_subtyped_count = 2, # H3N2, H1N1
    influenza_other_count = 2     # influenzaB, other
  )
}

# Helper function to get expected pathogen names
get_expected_pathogen_names <- function() {
  list(
    sarscov2_multiple = c("alpha", "delta", "omicron", "other"),
    influenza_subtyped = c("H3N2", "H1N1"),
    influenza_other = c("influenzaB", "other")
  )
}
