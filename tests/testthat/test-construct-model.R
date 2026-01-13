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
    expect_named(result, c("data", "validated_tsbl", "standata",
                           "pathogen_names", "dow_effect"))
    expect_s3_class(result, c(model_name, "EpiStrainDynamics.model"))

    # Test data structure
    expect_named(result$data, c("case_timeseries"))
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
    expect_equal(result$pathogen_names, expected_names$influenza_subtyped)
    expect_true("influenzaA_subtyped" %in% names(result$data))
  }
})

test_that("p_spline models include spline-specific parameters in standata", {
  check_package_data()

  models <- create_test_models()

  # Test that p-spline models have knots in standata while random walk models don't
  for (structure_type in c("single", "multiple", "subtyped")) {
    ps_model <- models[[paste0("ps_", structure_type)]]
    rw_model <- models[[paste0("rw_", structure_type)]]

    expect_true("knots" %in% names(ps_model$standata),
                info = paste("P-spline", structure_type, "should have knots in standata"))
    expect_false("knots" %in% names(rw_model$standata),
                 info = paste("Random walk", structure_type, "should not have knots in standata"))

    # Test spline parameters in standata
    expect_true("spline_degree" %in% names(ps_model$standata))
    expect_true("num_knots" %in% names(ps_model$standata))
    expect_true("X" %in% names(ps_model$standata))
  }
})

test_that("construct_model() handles day-of-week effects correctly", {
  check_package_data()

  method <- random_walk()
  pathogen <- single(
    data = sarscov2,
    case_timeseries = 'cases',
    time = 'date'
  )

  expected_length <- get_expected_data_lengths()$sarscov2_length

  # Test with dow_effect = TRUE
  result_dow <- construct_model(method, pathogen, dow_effect = TRUE)
  expect_equal(result_dow$standata$week_effect, 7L)
  expect_equal(result_dow$standata$DOW, ((1:expected_length - 1L) %% 7L) + 1L)
  expect_true(result_dow$dow_effect)

  # Test with dow_effect = FALSE (default)
  result_no_dow <- construct_model(method, pathogen, dow_effect = FALSE)
  expect_equal(result_no_dow$standata$week_effect, 1L)
  expect_equal(result_no_dow$standata$DOW, rep(1L, expected_length))
  expect_false(result_no_dow$dow_effect)
})

test_that("construct_model() validates inputs appropriately", {
  check_package_data()

  method <- random_walk()
  pathogen <- single(
    data = sarscov2,
    case_timeseries = 'cases',
    time = 'date'
  )

  # Test invalid method class
  invalid_method <- list(method = "random-walk")
  expect_error(construct_model(invalid_method, pathogen),
               "EpiStrainDynamics.method")

  # Test invalid pathogen structure class
  invalid_pathogen <- list(pathogen_structure = "single")
  expect_error(construct_model(method, invalid_pathogen),
               "EpiStrainDynamics.pathogen_structure")

  # Test dow_effect validation
  expect_error(construct_model(method, pathogen, dow_effect = "TRUE"),
               "must be a single logical value")

  # Test pathogen_noise validation
  expect_error(construct_model(method, pathogen, pathogen_noise = "TRUE"),
               "must be a single logical value")
})

test_that("smoothing parameters are correctly incorporated into standata", {
  check_package_data()

  method <- random_walk()
  pathogen <- multiple(
    data = sarscov2,
    case_timeseries = 'cases',
    time = 'date',
    component_pathogen_timeseries = c('alpha', 'delta', 'omicron', 'other')
  )

  # Test shared smoothing (default)
  result_shared <- construct_model(
    method, pathogen,
    smoothing_params = smoothing_structure("shared")
  )
  expect_equal(result_shared$standata$cov_structure, 0)
  expect_equal(result_shared$standata$tau_priors_provided, 1)

  # Test independent smoothing with custom priors
  result_indep <- construct_model(
    method, pathogen,
    smoothing_params = smoothing_structure("independent",
                                           tau_mean = c(0, 0.1, 0.3, 0),
                                           tau_sd = rep(1, 4))
  )
  expect_equal(result_indep$standata$cov_structure, 1)
  expect_equal(result_indep$standata$tau_priors_provided, 2)
  expect_equal(result_indep$standata$tau_mean, c(0, 0.1, 0.3, 0))
  expect_equal(result_indep$standata$tau_sd, rep(1, 4))

  # Test correlated smoothing
  result_corr <- construct_model(
    method, pathogen,
    smoothing_params = smoothing_structure("correlated")
  )
  expect_equal(result_corr$standata$cov_structure, 2)
})

test_that("dispersion parameters are correctly incorporated into standata", {
  check_package_data()

  method <- random_walk()
  pathogen <- single(
    data = sarscov2,
    case_timeseries = 'cases',
    time = 'date'
  )

  # Test default dispersion (no priors)
  result_default <- construct_model(
    method, pathogen,
    dispersion_params = dispersion_structure()
  )
  expect_equal(result_default$standata$phi_priors_provided, 1)

  # Test custom dispersion priors
  result_custom <- construct_model(
    method, pathogen,
    dispersion_params = dispersion_structure(phi_mean = 2.0, phi_sd = 0.5)
  )
  expect_equal(result_custom$standata$phi_priors_provided, 2)
  expect_equal(result_custom$standata$phi_mean, 2.0)
  expect_equal(result_custom$standata$phi_sd, 0.5)
})

test_that("pathogen_noise parameter is correctly incorporated into standata", {
  check_package_data()

  method <- random_walk()
  pathogen <- multiple(
    data = sarscov2,
    case_timeseries = 'cases',
    time = 'date',
    component_pathogen_timeseries = c('alpha', 'delta')
  )

  # Test pathogen_noise = FALSE (default)
  result_no_noise <- construct_model(method, pathogen, pathogen_noise = FALSE)
  expect_equal(result_no_noise$standata$noise_structure, 0)

  # Test pathogen_noise = TRUE
  result_with_noise <- construct_model(method, pathogen, pathogen_noise = TRUE)
  expect_equal(result_with_noise$standata$noise_structure, 1)
})

test_that("pathogen structure parameters are preserved correctly in standata", {
  check_package_data()

  models <- create_test_models()

  # Test multiple pathogen models have correct standata
  multiple_models <- list(
    rw_multiple = models$rw_multiple,
    ps_multiple = models$ps_multiple
  )

  for (model in multiple_models) {
    expect_true("cov_structure" %in% names(model$standata))
    expect_true("noise_structure" %in% names(model$standata))
    expect_true("num_path" %in% names(model$standata))
    expect_true("P" %in% names(model$standata))
  }

  # Test subtyped pathogen models have correct standata
  subtyped_models <- list(
    rw_subtyped = models$rw_subtyped,
    ps_subtyped = models$ps_subtyped
  )

  for (model in subtyped_models) {
    expect_true("cov_structure" %in% names(model$standata))
    expect_true("noise_structure" %in% names(model$standata))
    expect_true("num_path" %in% names(model$standata))
    expect_true("P1" %in% names(model$standata))
    expect_true("P2" %in% names(model$standata))

    # Verify data is also in model$data
    expect_true("component_pathogens" %in% names(model$data))
    expect_true("influenzaA_subtyped" %in% names(model$data))
  }
})

test_that("models with various parameter combinations work correctly", {
  check_package_data()

  models <- create_test_models()

  # Test model with independent smoothing
  expect_equal(models$rw_multiple_indep_smooth$standata$cov_structure, 1)
  expect_equal(models$rw_multiple_indep_smooth$standata$tau_mean, c(0, 0.1, 0.3, 0))

  # Test model with correlated smoothing
  expect_equal(models$rw_multiple_corr_smooth$standata$cov_structure, 2)

  # Test model with custom dispersion
  expect_equal(models$rw_multiple_custom_disp$standata$phi_mean, 2.0)
  expect_equal(models$rw_multiple_custom_disp$standata$phi_sd, 0.5)

  # Test model with pathogen noise
  expect_equal(models$rw_multiple_noise$standata$noise_structure, 1)

  # Test model with dow effect
  expect_equal(models$rw_single_dow$standata$week_effect, 7L)
  expect_true(models$rw_single_dow$dow_effect)

  # Test full-featured models
  expect_equal(models$rw_multiple_full$standata$cov_structure, 1)
  expect_equal(models$rw_multiple_full$standata$noise_structure, 1)
  expect_equal(models$rw_multiple_full$standata$week_effect, 7L)
  expect_true(models$rw_multiple_full$dow_effect)
})

# Tests for get_model_type() - kept as is since they test pure logic
test_that("get_model_type() returns correct model types", {
  expect_equal(EpiStrainDynamics:::get_model_type("random-walk", "single"), "rw_single")
  expect_equal(EpiStrainDynamics:::get_model_type("random-walk", "multiple"), "rw_multiple")
  expect_equal(EpiStrainDynamics:::get_model_type("random-walk", "subtyped"), "rw_subtyped")
  expect_equal(EpiStrainDynamics:::get_model_type("p-spline", "single"), "ps_single")
  expect_equal(EpiStrainDynamics:::get_model_type("p-spline", "multiple"), "ps_multiple")
  expect_equal(EpiStrainDynamics:::get_model_type("p-spline", "subtyped"), "ps_subtyped")
})

test_that("get_model_type() validates inputs and provides helpful errors", {
  # Test invalid methods
  expect_error(EpiStrainDynamics:::get_model_type("invalid-method", "single"),
               "Unknown `method`.*invalid-method")
  expect_error(EpiStrainDynamics:::get_model_type("random_walk", "single"),
               "Unknown `method`.*random_walk")

  # Test invalid structures
  expect_error(EpiStrainDynamics:::get_model_type("random-walk", "invalid-structure"),
               "Unknown `pathogen_type`.*invalid-structure")
  expect_error(EpiStrainDynamics:::get_model_type("p-spline", ""),
               "Unknown `pathogen_type`")
})

# Tests for get_knots() - kept similar but updated for consistency
test_that("get_knots() calculates knots correctly", {
  X <- 1:10
  knots <- EpiStrainDynamics:::get_knots(X, days_per_knot = 2, spline_degree = 1)

  expect_type(knots, "double")
  expect_true(length(knots) > 0)
  expect_true(all(is.finite(knots)))
  expect_true(min(knots) < min(X))
  expect_true(max(knots) > max(X))
})

test_that("get_knots() handles parameter variations correctly", {
  X <- 1:20

  # More knots for smaller days_per_knot
  knots_2 <- EpiStrainDynamics:::get_knots(X, days_per_knot = 2)
  knots_5 <- EpiStrainDynamics:::get_knots(X, days_per_knot = 5)
  expect_true(length(knots_2) > length(knots_5))

  # Wider range for higher spline degree
  knots_deg1 <- EpiStrainDynamics:::get_knots(X, spline_degree = 1)
  knots_deg4 <- EpiStrainDynamics:::get_knots(X, spline_degree = 4)
  expect_true(min(knots_deg4) < min(knots_deg1))
  expect_true(max(knots_deg4) > max(knots_deg1))
})

test_that("get_knots() validates inputs properly", {
  X <- 1:10

  # Test invalid X
  expect_error(EpiStrainDynamics:::get_knots(character(5)), "must be a non-empty numeric vector")
  expect_error(EpiStrainDynamics:::get_knots(numeric(0)), "must be a non-empty numeric vector")

  # Test invalid parameters
  expect_error(EpiStrainDynamics:::get_knots(X, days_per_knot = 0), "must be a positive")
  expect_error(EpiStrainDynamics:::get_knots(X, days_per_knot = 2.5), "must be a whole number")
  expect_error(EpiStrainDynamics:::get_knots(X, spline_degree = 0), "must be a positive")
  expect_error(EpiStrainDynamics:::get_knots(X, spline_degree = 3.7), "must be a whole number")
})

test_that("get_cov_structure() returns correct numeric codes", {
  expect_equal(EpiStrainDynamics:::get_cov_structure("shared"), 0)
  expect_equal(EpiStrainDynamics:::get_cov_structure("independent"), 1)
  expect_equal(EpiStrainDynamics:::get_cov_structure("correlated"), 2)

  # Test error handling
  expect_error(EpiStrainDynamics:::get_cov_structure("invalid"), "Invalid option provided")
})

# Integration tests using helper functions
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

    # Test method-specific features in standata
    if (startsWith(test_case$name, "ps_")) {
      expect_true("knots" %in% names(model$standata))
      expect_true("spline_degree" %in% names(model$standata))
      expect_true("num_knots" %in% names(model$standata))
    } else {
      expect_false("knots" %in% names(model$standata))
    }

    # Test required standata components
    expect_true(all(c("num_data", "Y", "week_effect", "DOW") %in% names(model$standata)))
  }
})

test_that("Custom model construction with specific parameters", {
  check_package_data()

  # Test custom p-spline with specific parameters
  method <- p_spline(spline_degree = 4, days_per_knot = 2)
  pathogen <- multiple(
    data = sarscov2,
    case_timeseries = 'cases',
    time = 'date',
    component_pathogen_timeseries = c('alpha', 'delta', 'omicron')
  )

  result <- construct_model(
    method,
    pathogen,
    smoothing_params = smoothing_structure("correlated"),
    pathogen_noise = TRUE,
    dow_effect = TRUE
  )

  expect_s3_class(result, c("ps_multiple", "EpiStrainDynamics.model"))
  expect_equal(result$standata$spline_degree, 4L)
  expect_equal(result$standata$week_effect, 7L)
  expect_equal(result$standata$cov_structure, 2)  # correlated
  expect_equal(result$standata$noise_structure, 1)  # pathogen noise = TRUE
  expect_equal(result$pathogen_names, c('alpha', 'delta', 'omicron'))
  expect_true(result$dow_effect)
})
