# Tests for random_walk() function
test_that("random_walk() returns correct structure and class", {
  result <- random_walk()

  # Test that result is a list
  expect_type(result, "list")

  # Test that it has the correct class
  expect_s3_class(result, "EpiStrainDynamics.method")

  # Test that it contains the correct method name
  expect_equal(result$method, "random-walk")

  # Test that the list has only one element
  expect_length(result, 1)

  # Test that method element exists and is character
  expect_true("method" %in% names(result))
  expect_type(result$method, "character")
})

test_that("random_walk() returns consistent results", {
  result1 <- random_walk()
  result2 <- random_walk()

  # Should return identical results every time (no randomness in function itself)
  expect_identical(result1, result2)
})

# Tests for p_spline() function
test_that("p_spline() with default parameters returns correct structure", {
  result <- p_spline()

  # Test that result is a list
  expect_type(result, "list")

  # Test that it has the correct class
  expect_s3_class(result, "EpiStrainDynamics.method")

  # Test that it contains the correct method name
  expect_equal(result$method, "p-spline")

  # Test that it has model_params
  expect_true("model_params" %in% names(result))
  expect_type(result$model_params, "list")

  # Test default parameter values
  expect_equal(result$model_params$spline_degree, 3)
  expect_equal(result$model_params$days_per_knot, 3)

  # Test structure has exactly 2 top-level elements
  expect_length(result, 2)
  expect_length(result$model_params, 2)
})

test_that("p_spline() with custom parameters works correctly", {
  result <- p_spline(spline_degree = 2, days_per_knot = 5)

  # Test custom parameter values are set correctly
  expect_equal(result$model_params$spline_degree, 2)
  expect_equal(result$model_params$days_per_knot, 5)

  # Test that method is still correct
  expect_equal(result$method, "p-spline")

  # Test class is still correct
  expect_s3_class(result, "EpiStrainDynamics.method")
})

test_that("p_spline() parameter types are preserved", {
  result <- p_spline(spline_degree = 1L, days_per_knot = 10L)

  # Test that numeric inputs are preserved as numeric
  expect_type(result$model_params$spline_degree, "integer")
  expect_type(result$model_params$days_per_knot, "integer")
})

test_that("p_spline() handles edge case values", {
  # Test with minimum reasonable values
  result_min <- p_spline(spline_degree = 1, days_per_knot = 1)
  expect_equal(result_min$model_params$spline_degree, 1)
  expect_equal(result_min$model_params$days_per_knot, 1)

  # Test with larger values
  result_large <- p_spline(spline_degree = 10, days_per_knot = 30)
  expect_equal(result_large$model_params$spline_degree, 10)
  expect_equal(result_large$model_params$days_per_knot, 30)
})

test_that("Both functions return objects with same class", {
  rw_result <- random_walk()
  ps_result <- p_spline()

  # Both should have the same class
  expect_identical(class(rw_result), class(ps_result))
  expect_s3_class(rw_result, "EpiStrainDynamics.method")
  expect_s3_class(ps_result, "EpiStrainDynamics.method")
})

test_that("Function outputs have expected names", {
  rw_result <- random_walk()
  ps_result <- p_spline()

  # random_walk should only have 'method'
  expect_named(rw_result, "method")

  # p_spline should have 'method' and 'model_params'
  expect_named(ps_result, c("method", "model_params"))
  expect_named(ps_result$model_params, c("spline_degree", "days_per_knot"))
})

# Additional validation tests
test_that("p_spline() accepts decimal values", {
  result <- p_spline(spline_degree = 2.5, days_per_knot = 3.7)
  expect_equal(result$model_params$spline_degree, 2.5)
  expect_equal(result$model_params$days_per_knot, 3.7)
})
