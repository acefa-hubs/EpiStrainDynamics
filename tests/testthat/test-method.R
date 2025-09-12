# Tests for method constructor functions - simplified version
#' @srrstats {G3.0} comparisons made between appropriate values

# Helper function for method validation tests (add to helper script)
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

# Tests for random_walk() function
test_that("random_walk() returns correct structure and class", {
  result <- random_walk()
  test_method_structure(result, "random-walk")
  expect_named(result, "method")
})

# Tests for p_spline() function
test_that("p_spline() with default parameters returns correct structure", {
  result <- p_spline()
  expected_params <- list(spline_degree = 3, days_per_knot = 3)

  test_method_structure(result, "p-spline", expected_params)
  expect_named(result, c("method", "model_params"))
  expect_length(result, 2)
  expect_length(result$model_params, 2)
})

test_that("p_spline() with custom parameters works correctly", {
  result <- p_spline(spline_degree = 2, days_per_knot = 5)
  expected_params <- list(spline_degree = 2, days_per_knot = 5)

  test_method_structure(result, "p-spline", expected_params)
})

test_that("p_spline() handles parameter types and edge cases correctly", {
  # Test type preservation
  result_int <- p_spline(spline_degree = 1L, days_per_knot = 10L)
  expect_type(result_int$model_params$spline_degree, "integer")
  expect_type(result_int$model_params$days_per_knot, "integer")

  # Test conversion of whole numbers
  result_converted <- p_spline(spline_degree = 3.0, days_per_knot = 5.0)
  expect_type(result_converted$model_params$spline_degree, "integer")
  expect_type(result_converted$model_params$days_per_knot, "integer")
  expect_equal(result_converted$model_params$spline_degree, 3L)
  expect_equal(result_converted$model_params$days_per_knot, 5L)

  # Test edge case values
  result_min <- p_spline(spline_degree = 1, days_per_knot = 1)
  expected_min <- list(spline_degree = 1, days_per_knot = 1)
  test_method_structure(result_min, "p-spline", expected_min)

  result_large <- p_spline(spline_degree = 10, days_per_knot = 30)
  expected_large <- list(spline_degree = 10, days_per_knot = 30)
  test_method_structure(result_large, "p-spline", expected_large)
})

test_that("Both method constructors return compatible objects", {
  rw_result <- random_walk()
  ps_result <- p_spline()

  # Both should have the same class
  expect_identical(class(rw_result), class(ps_result))
  expect_s3_class(rw_result, "EpiStrainDynamics.method")
  expect_s3_class(ps_result, "EpiStrainDynamics.method")
})

# Validation tests - these stay mostly the same as they test specific validation logic
test_that("p_spline() validates numeric inputs correctly", {
  # Non-numeric inputs
  expect_error(p_spline(spline_degree = "3", days_per_knot = 3),
               "must be numeric")
  expect_error(p_spline(spline_degree = 3, days_per_knot = "5"),
               "must be numeric")

  # Invalid special values
  expect_error(p_spline(spline_degree = NA, days_per_knot = 3),
               "must be numeric")
  expect_error(p_spline(spline_degree = 3, days_per_knot = Inf),
               "must be a finite number")
  expect_error(p_spline(spline_degree = NULL, days_per_knot = 3),
               "must be numeric")

  # Vector inputs (length > 1)
  expect_error(p_spline(spline_degree = c(2, 3), days_per_knot = 3),
               "must be a single value")
  expect_error(p_spline(spline_degree = 3, days_per_knot = c(3, 4)),
               "must be a single value")
})

test_that("p_spline() validates positive whole numbers correctly", {
  # Test decimal values are rejected
  expect_error(p_spline(spline_degree = 2.5, days_per_knot = 3),
               "must be a whole number")
  expect_error(p_spline(spline_degree = 3, days_per_knot = 3.7),
               "must be a whole number")
  expect_error(p_spline(spline_degree = 2.1, days_per_knot = 4.9),
               "must be a whole number")

  # Test zero and negative values are rejected
  expect_error(p_spline(spline_degree = 0, days_per_knot = 3),
               "must be a positive number")
  expect_error(p_spline(spline_degree = 3, days_per_knot = 0),
               "must be a positive number")
  expect_error(p_spline(spline_degree = -1, days_per_knot = 3),
               "must be a positive number")
  expect_error(p_spline(spline_degree = 3, days_per_knot = -5),
               "must be a positive number")
})
