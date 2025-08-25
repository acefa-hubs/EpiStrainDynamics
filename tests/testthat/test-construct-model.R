# Tests for construct_model() with real functions
test_that("construct_model() works with random_walk() and single()", {

  method <- random_walk()
  pathogen <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "TestPathogen"
  )

  result <- construct_model(method, pathogen, dow_effect = FALSE)

  # Test structure
  expect_type(result, "list")
  expect_named(result, c("data", "model_params", "pathogen_names", "dow_effect"))

  # Test classes
  expect_s3_class(result, "rw_single")
  expect_s3_class(result, "EpiStrainDynamics.model")

  # Test data structure
  expect_named(result$data, c("time_seq", "case_timeseries", "time"))
  expect_equal(result$data$time_seq, 1:nrow(sarscov2))
  expect_equal(result$data$case_timeseries, sarscov2$cases)
  expect_equal(result$data$time, sarscov2$date)
  expect_equal(result$pathogen_names, "TestPathogen")
})

test_that("construct_model() works with p_spline() and multiple()", {

  method <- p_spline(spline_degree = 2, days_per_knot = 4)
  pathogen <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta,
      omicron = sarscov2$omicron
    ),
    smoothing_structure = "independent",
    observation_noise = "observation_noise_only"
  )

  result <- construct_model(method, pathogen, dow_effect = TRUE)

  # Test classes
  expect_s3_class(result, "ps_multiple")
  expect_s3_class(result, "EpiStrainDynamics.model")

  # Test method parameters are included
  expect_equal(result$model_params$spline_degree, 2L)
  expect_equal(result$model_params$days_per_knot, 4L)

  # Test that knots are calculated
  expect_true("knots" %in% names(result$model_params))
  expect_type(result$model_params$knots, "double")

  # Test pathogen structure parameters
  expect_equal(result$model_params$cov_structure, 1)  # independent
  expect_equal(result$model_params$noise_structure, 0)  # observation_noise_only

  # Test pathogen names
  expect_equal(result$pathogen_names, c("alpha", "delta", "omicron"))

  # Test day-of-week effect
  expect_equal(result$model_params$week_effect, 7L)
  expect_true(result$dow_effect)
})

test_that("construct_model() works with subtyped pathogen structure", {

  method <- random_walk()
  pathogen <- subtyped(
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
    ),
    smoothing_structure = "correlated",
    observation_noise = "pathogen_specific_noise"
  )

  result <- construct_model(method, pathogen, dow_effect = FALSE)

  # Test classes
  expect_s3_class(result, "rw_subtyped")
  expect_s3_class(result, "EpiStrainDynamics.model")

  # Test pathogen structure parameters
  expect_equal(result$model_params$cov_structure, 2)  # correlated
  expect_equal(result$model_params$noise_structure, 1)  # pathogen_specific_noise

  # Test pathogen names (from subtyped and other pathogens)
  expect_equal(result$pathogen_names, c("H3N2", "H1N1", "influenzaB", "other"))

  # Test data structure includes subtyped data
  expect_true("component_pathogens" %in% names(result$data))
  expect_true("influenzaA_subtyped" %in% names(result$data))
})

test_that("construct_model() handles day-of-week effects with real functions", {

  method <- random_walk()
  pathogen <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "TestPathogen"
  )

  # Test with dow_effect = TRUE
  result_dow <- construct_model(method, pathogen, dow_effect = TRUE)
  expect_equal(result_dow$model_params$week_effect, 7L)
  expect_equal(result_dow$model_params$DOW, ((1:nrow(sarscov2) - 1L) %% 7L) + 1L)
  expect_true(result_dow$dow_effect)

  # Test with dow_effect = FALSE
  result_no_dow <- construct_model(method, pathogen, dow_effect = FALSE)
  expect_equal(result_no_dow$model_params$week_effect, 1L)
  expect_equal(result_no_dow$model_params$DOW, rep(1L, nrow(sarscov2)))
  expect_false(result_no_dow$dow_effect)
})

test_that("construct_model() validates inputs from real functions", {

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

test_that("construct_model() preserves all pathogen structure data", {

  method <- p_spline(spline_degree = 3, days_per_knot = 2)
  pathogen <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta,
      omicron = sarscov2$omicron
    ),
    smoothing_structure = "shared",
    observation_noise = "pathogen_specific_noise"
  )

  result <- construct_model(method, pathogen, dow_effect = FALSE)

  # Test that all original data is preserved
  expect_equal(result$data$case_timeseries, sarscov2$cases)
  expect_equal(result$data$time, sarscov2$date)
  expect_true("component_pathogens" %in% names(result$data))

  # Test that pathogen model params are preserved
  expect_equal(result$model_params$cov_structure, 0)  # shared
  expect_equal(result$model_params$noise_structure, 1)  # pathogen_specific_noise

  # Test pathogen names
  expect_equal(result$pathogen_names, c("alpha", "delta", "omicron"))
})

test_that("p_spline validation is applied in real integration", {

  pathogen <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "TestPathogen"
  )

  # Test that p_spline validation works
  expect_error(p_spline(spline_degree = 2.5), "must be a whole number")
  expect_error(p_spline(days_per_knot = 0), "must be a positive number")

  # Test that valid p_spline works in construct_model
  method <- p_spline(spline_degree = 4, days_per_knot = 2)
  result <- construct_model(method, pathogen, dow_effect = FALSE)

  expect_equal(result$model_params$spline_degree, 4L)
  expect_equal(result$model_params$days_per_knot, 2L)
})

# Tests for get_model_type()
test_that("get_model_type() returns correct model types", {
  expect_equal(get_model_type("random-walk", "single"), "rw_single")
  expect_equal(get_model_type("random-walk", "multiple"), "rw_multiple")
  expect_equal(get_model_type("random-walk", "subtyped"), "rw_subtyped")
  expect_equal(get_model_type("p-spline", "single"), "ps_single")
  expect_equal(get_model_type("p-spline", "multiple"), "ps_multiple")
  expect_equal(get_model_type("p-spline", "subtyped"), "ps_subtyped")
})

test_that("get_model_type() validates method names", {
  expect_error(get_model_type("invalid-method", "single"),
               "Unknown method: invalid-method")
  expect_error(get_model_type("random_walk", "single"),
               "Unknown method: random_walk")
  expect_error(get_model_type("", "single"),
               "Unknown method:")
})

test_that("get_model_type() validates pathogen structure types", {
  expect_error(get_model_type("random-walk", "invalid-structure"),
               "Unknown pathogen structure: invalid-structure")
  expect_error(get_model_type("p-spline", ""),
               "Unknown pathogen structure:")
})

test_that("get_model_type() provides helpful error messages", {
  expect_error(get_model_type("invalid", "single"),
               "Valid methods are: random-walk, p-spline")
  expect_error(get_model_type("random-walk", "invalid"),
               "Valid structures are: single, multiple, subtyped")
})

# Tests for get_knots()
test_that("get_knots() calculates knots correctly", {
  X <- 1:10
  knots <- get_knots(X, days_per_knot = 2, spline_degree = 1)

  expect_type(knots, "double")
  expect_true(length(knots) > 0)
  expect_true(all(is.finite(knots)))

  # Test that knots extend beyond the data range
  expect_true(min(knots) < min(X))
  expect_true(max(knots) > max(X))
})

test_that("get_knots() uses default parameters correctly", {
  X <- 1:15
  knots_default <- get_knots(X)
  knots_explicit <- get_knots(X, days_per_knot = 3, spline_degree = 3)

  expect_equal(knots_default, knots_explicit)
})

test_that("get_knots() handles different parameter values", {
  X <- 1:20

  # Test different days_per_knot
  knots_2 <- get_knots(X, days_per_knot = 2)
  knots_5 <- get_knots(X, days_per_knot = 5)
  expect_true(length(knots_2) > length(knots_5))

  # Test different spline_degree
  knots_deg1 <- get_knots(X, spline_degree = 1)
  knots_deg4 <- get_knots(X, spline_degree = 4)
  expect_true(min(knots_deg4) < min(knots_deg1))
  expect_true(max(knots_deg4) > max(knots_deg1))
})

test_that("get_knots() validates input parameters", {
  X <- 1:10

  # Test invalid X
  expect_error(get_knots(character(5)), "must be a non-empty numeric vector")
  expect_error(get_knots(numeric(0)), "must be a non-empty numeric vector")
  expect_error(get_knots(NULL), "must be a non-empty numeric vector")

  # Test invalid days_per_knot
  expect_error(get_knots(X, days_per_knot = 0), "must be a positive number")
  expect_error(get_knots(X, days_per_knot = -1), "must be a positive number")
  expect_error(get_knots(X, days_per_knot = 2.5), "must be a whole number")
  expect_error(get_knots(X, days_per_knot = "3"), "must be numeric")

  # Test invalid spline_degree
  expect_error(get_knots(X, spline_degree = 0), "must be a positive number")
  expect_error(get_knots(X, spline_degree = -2), "must be a positive number")
  expect_error(get_knots(X, spline_degree = 3.7), "must be a whole number")
})

test_that("get_knots() handles edge cases", {
  # Test single value
  expect_silent(get_knots(5))

  # Test with NA values
  X_with_na <- c(1, 2, NA, 4, 5)
  result <- get_knots(X_with_na)
  expect_type(result, "double")
  expect_true(all(is.finite(result)))

  # Test all NA values
  expect_error(get_knots(c(NA, NA, NA)), "must be a non-empty numeric vector")
})

test_that("get_knots() handles non-integer sequences", {
  X <- c(1.5, 2.7, 3.1, 4.8, 5.9)
  result <- get_knots(X, days_per_knot = 1, spline_degree = 1)

  expect_type(result, "double")
  expect_true(length(result) > 0)
  expect_true(min(result) < min(X))
  expect_true(max(result) > max(X))
})

# Integration tests with real functions
test_that("Full integration: random walk with single pathogen", {

  method <- random_walk()
  pathogen <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "SARS-CoV-2"
  )

  result <- construct_model(method, pathogen, dow_effect = TRUE)

  expect_s3_class(result, c("rw_single", "EpiStrainDynamics.model", "list"))
  expect_equal(length(result$data$time_seq), nrow(sarscov2))
  expect_equal(result$model_params$week_effect, 7L)
  expect_false("knots" %in% names(result$model_params))
  expect_equal(result$pathogen_names, "SARS-CoV-2")
})

test_that("Full integration: p-spline with multiple pathogens", {

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

  result <- construct_model(method, pathogen, dow_effect = FALSE)

  expect_s3_class(result, c("ps_multiple", "EpiStrainDynamics.model", "list"))
  expect_equal(result$model_params$spline_degree, 4L)
  expect_equal(result$model_params$days_per_knot, 2L)
  expect_true("knots" %in% names(result$model_params))
  expect_equal(result$model_params$week_effect, 1L)
  expect_equal(result$model_params$cov_structure, 2)  # correlated
  expect_equal(result$model_params$noise_structure, 1)  # pathogen_specific_noise
  expect_equal(result$pathogen_names, c("alpha", "delta", "omicron"))
})

test_that("Full integration: p-spline with subtyped pathogens", {

  method <- p_spline(spline_degree = 2, days_per_knot = 3)
  pathogen <- subtyped(
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
    ),
    smoothing_structure = "independent",
    observation_noise = "observation_noise_only"
  )

  result <- construct_model(method, pathogen, dow_effect = TRUE)

  expect_s3_class(result, c("ps_subtyped", "EpiStrainDynamics.model", "list"))
  expect_equal(result$model_params$spline_degree, 2L)
  expect_equal(result$model_params$days_per_knot, 3L)
  expect_true("knots" %in% names(result$model_params))
  expect_equal(result$model_params$week_effect, 7L)
  expect_equal(result$model_params$cov_structure, 1)  # independent
  expect_equal(result$model_params$noise_structure, 0)  # observation_noise_only
  expect_equal(result$pathogen_names, c("H3N2", "H1N1", "influenzaB", "other"))

  # Test subtyped-specific data structure
  expect_true("influenzaA_subtyped" %in% names(result$data))
  expect_equal(nrow(result$data$influenzaA_subtyped), 2)  # 2 subtypes
  expect_equal(ncol(result$data$influenzaA_subtyped), nrow(influenza))  # 10 time points
})
