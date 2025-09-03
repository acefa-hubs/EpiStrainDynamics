# Revised tests using helper functions for better maintainability

# Tests for construct_model() using standardized test models
test_that("construct_model() works with all standard model configurations", {
  check_package_data()

  models <- create_test_models()
  expected_lengths <- get_expected_data_lengths()
  expected_names <- get_expected_pathogen_names()

  # Test single pathogen models
  for (method_type in c("rw", "ps")) {
    model_name <- paste0(method_type, "_single")
    result <- models[[model_name]]

    # Test structure and classes
    expect_type(result, "list")
    expect_named(result, c("data", "model_params", "pathogen_names", "dow_effect"))
    expect_s3_class(result, c(model_name, "EpiStrainDynamics.model"))

    # Test data structure
    expect_named(result$data, c("time_seq", "case_timeseries", "time"))
    expect_equal(result$data$time_seq, 1:expected_lengths$sarscov2_length)
    expect_equal(length(result$data$case_timeseries), expected_lengths$sarscov2_length)
  }

  # Test multiple pathogen models
  for (method_type in c("rw", "ps")) {
    model_name <- paste0(method_type, "_multiple")
    result <- models[[model_name]]

    expect_s3_class(result, c(model_name, "EpiStrainDynamics.model"))
    expect_equal(result$pathogen_names, expected_names$sarscov2_multiple)
    expect_true("component_pathogens" %in% names(result$data))
  }

  # Test subtyped pathogen models
  for (method_type in c("rw", "ps")) {
    model_name <- paste0(method_type, "_subtyped")
    result <- models[[model_name]]

    expect_s3_class(result, c(model_name, "EpiStrainDynamics.model"))
    expect_equal(result$pathogen_names, c(expected_names$influenza_subtyped, expected_names$influenza_other))
    expect_true("influenzaA_subtyped" %in% names(result$data))
  }
})

test_that("p_spline models include spline-specific parameters", {
  check_package_data()

  models <- create_test_models()

  # Test that p-spline models have knots while random walk models don't
  for (structure_type in c("single", "multiple", "subtyped")) {
    ps_model <- models[[paste0("ps_", structure_type)]]
    rw_model <- models[[paste0("rw_", structure_type)]]

    expect_true("knots" %in% names(ps_model$model_params),
                info = paste("P-spline", structure_type, "should have knots"))
    expect_false("knots" %in% names(rw_model$model_params),
                 info = paste("Random walk", structure_type, "should not have knots"))

    # Test spline parameters
    expect_true("spline_degree" %in% names(ps_model$model_params))
    expect_true("days_per_knot" %in% names(ps_model$model_params))
  }
})

test_that("construct_model() handles day-of-week effects correctly", {
  check_package_data()

  method <- random_walk()
  pathogen <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "TestPathogen"
  )

  expected_length <- get_expected_data_lengths()$sarscov2_length

  # Test with dow_effect = TRUE
  result_dow <- construct_model(method, pathogen, dow_effect = TRUE)
  expect_equal(result_dow$model_params$week_effect, 7L)
  expect_equal(result_dow$model_params$DOW, ((1:expected_length - 1L) %% 7L) + 1L)
  expect_true(result_dow$dow_effect)

  # Test with dow_effect = FALSE
  result_no_dow <- construct_model(method, pathogen, dow_effect = FALSE)
  expect_equal(result_no_dow$model_params$week_effect, 1L)
  expect_equal(result_no_dow$model_params$DOW, rep(1L, expected_length))
  expect_false(result_no_dow$dow_effect)
})

test_that("construct_model() validates inputs appropriately", {
  check_package_data()

  method <- random_walk()
  pathogen <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "TestPathogen"
  )

  # Test invalid method class
  invalid_method <- list(method = "random-walk")
  expect_error(construct_model(invalid_method, pathogen),
               "Input must be of class EpiStrainDynamics.method")

  # Test invalid pathogen structure class
  invalid_pathogen <- list(pathogen_structure = "single")
  expect_error(construct_model(method, invalid_pathogen),
               "Input must be of class EpiStrainDynamics.pathogen_structure")

  # Test dow_effect validation
  expect_error(construct_model(method, pathogen, dow_effect = "TRUE"),
               "must be a single logical value")
})

test_that("pathogen structure parameters are preserved correctly", {
  check_package_data()

  models <- create_test_models()

  # Test multiple pathogen smoothing and noise structures
  multiple_models <- list(
    rw_multiple = models$rw_multiple,
    ps_multiple = models$ps_multiple
  )

  for (model in multiple_models) {
    expect_true("cov_structure" %in% names(model$model_params))
    expect_true("noise_structure" %in% names(model$model_params))
  }

  # Test subtyped pathogen structures
  subtyped_models <- list(
    rw_subtyped = models$rw_subtyped,
    ps_subtyped = models$ps_subtyped
  )

  for (model in subtyped_models) {
    expect_true("cov_structure" %in% names(model$model_params))
    expect_true("noise_structure" %in% names(model$model_params))
    expect_true("component_pathogens" %in% names(model$data))
    expect_true("influenzaA_subtyped" %in% names(model$data))
  }
})

# Tests for get_model_type() - kept as is since they test pure logic
test_that("get_model_type() returns correct model types", {
  expect_equal(get_model_type("random-walk", "single"), "rw_single")
  expect_equal(get_model_type("random-walk", "multiple"), "rw_multiple")
  expect_equal(get_model_type("random-walk", "subtyped"), "rw_subtyped")
  expect_equal(get_model_type("p-spline", "single"), "ps_single")
  expect_equal(get_model_type("p-spline", "multiple"), "ps_multiple")
  expect_equal(get_model_type("p-spline", "subtyped"), "ps_subtyped")
})

test_that("get_model_type() validates inputs and provides helpful errors", {
  # Test invalid methods
  expect_error(get_model_type("invalid-method", "single"),
               "Unknown method: invalid-method.*Valid methods are: random-walk, p-spline")
  expect_error(get_model_type("random_walk", "single"),
               "Unknown method: random_walk.*Valid methods are: random-walk, p-spline")

  # Test invalid structures
  expect_error(get_model_type("random-walk", "invalid-structure"),
               "Unknown pathogen structure: invalid-structure.*Valid structures are: single, multiple, subtyped")
  expect_error(get_model_type("p-spline", ""),
               "Unknown pathogen structure:.*Valid structures are: single, multiple, subtyped")
})

# Tests for get_knots() - simplified validation tests
test_that("get_knots() calculates knots correctly", {
  X <- 1:10
  knots <- get_knots(X, days_per_knot = 2, spline_degree = 1)

  expect_type(knots, "double")
  expect_true(length(knots) > 0)
  expect_true(all(is.finite(knots)))
  expect_true(min(knots) < min(X))
  expect_true(max(knots) > max(X))
})

test_that("get_knots() handles parameter variations correctly", {
  X <- 1:20

  # More knots for smaller days_per_knot
  knots_2 <- get_knots(X, days_per_knot = 2)
  knots_5 <- get_knots(X, days_per_knot = 5)
  expect_true(length(knots_2) > length(knots_5))

  # Wider range for higher spline degree
  knots_deg1 <- get_knots(X, spline_degree = 1)
  knots_deg4 <- get_knots(X, spline_degree = 4)
  expect_true(min(knots_deg4) < min(knots_deg1))
  expect_true(max(knots_deg4) > max(knots_deg1))
})

test_that("get_knots() validates inputs properly", {
  X <- 1:10

  # Test invalid X
  expect_error(get_knots(character(5)), "must be a non-empty numeric vector")
  expect_error(get_knots(numeric(0)), "must be a non-empty numeric vector")

  # Test invalid parameters
  expect_error(get_knots(X, days_per_knot = 0), "must be a positive number")
  expect_error(get_knots(X, days_per_knot = 2.5), "must be a whole number")
  expect_error(get_knots(X, spline_degree = 0), "must be a positive number")
  expect_error(get_knots(X, spline_degree = 3.7), "must be a whole number")
})

# Simplified integration tests using helper functions
test_that("Full integration: all model combinations work correctly", {
  check_package_data()

  models <- create_test_models()
  expected_lengths <- get_expected_data_lengths()
  expected_names <- get_expected_pathogen_names()

  # Test each model type
  test_cases <- list(
    list(name = "rw_single", pathogen_count = 1, dataset = "sarscov2"),
    list(name = "ps_single", pathogen_count = 1, dataset = "sarscov2"),
    list(name = "rw_multiple", pathogen_count = expected_lengths$sarscov2_pathogen_count, dataset = "sarscov2"),
    list(name = "ps_multiple", pathogen_count = expected_lengths$sarscov2_pathogen_count, dataset = "sarscov2"),
    list(name = "rw_subtyped", pathogen_count = expected_lengths$influenza_subtyped_count + expected_lengths$influenza_other_count, dataset = "influenza"),
    list(name = "ps_subtyped", pathogen_count = expected_lengths$influenza_subtyped_count + expected_lengths$influenza_other_count, dataset = "influenza")
  )

  for (test_case in test_cases) {
    model <- models[[test_case$name]]

    # Test basic structure
    expect_s3_class(model, c(test_case$name, "EpiStrainDynamics.model"))
    expect_length(model$pathogen_names, test_case$pathogen_count)

    # Test dataset-specific length
    expected_length <- if (test_case$dataset == "sarscov2") {
      expected_lengths$sarscov2_length
    } else {
      expected_lengths$influenza_length
    }
    expect_equal(length(model$data$time_seq), expected_length)

    # Test method-specific features
    if (startsWith(test_case$name, "ps_")) {
      expect_true("knots" %in% names(model$model_params))
      expect_true("spline_degree" %in% names(model$model_params))
    } else {
      expect_false("knots" %in% names(model$model_params))
    }
  }
})

test_that("Custom model construction with specific parameters", {
  check_package_data()

  # Test custom p-spline with specific parameters
  method <- p_spline(spline_degree = 4, days_per_knot = 2)
  pathogen <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta,
      omicron = sarscov2$omicron
    ),
    smoothing_structure = "correlated",
    observation_noise = "pathogen_specific_noise"
  )

  result <- construct_model(method, pathogen, dow_effect = TRUE)

  expect_s3_class(result, c("ps_multiple", "EpiStrainDynamics.model"))
  expect_equal(result$model_params$spline_degree, 4L)
  expect_equal(result$model_params$days_per_knot, 2L)
  expect_equal(result$model_params$week_effect, 7L)
  expect_equal(result$model_params$cov_structure, 2)  # correlated
  expect_equal(result$model_params$noise_structure, 1)  # pathogen_specific_noise
  expect_equal(result$pathogen_names, c('alpha', 'delta', 'omicron'))
})
