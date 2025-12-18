#' @srrstats {G5.2, G5.2b} error and warning behaviour explicitly tested,
#'   including conditions that trigger these messages and comparing with
#'   expected results
#' @srrstats {G3.0} comparisons made between appropriate values
#'
# Tests for validate_positive_whole_number()
test_that("validate_positive_whole_number() passes valid inputs", {
  # Should not throw errors for valid positive whole numbers
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(1, "test_param"))
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(5L, "test_param"))
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(100, "test_param"))
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(3.0, "test_param"))  # Whole number as double
})

test_that("validate_positive_whole_number() returns invisible NULL", {
  result <- EpiStrainDynamics:::validate_positive_whole_number(5, "test_param")
  expect_null(result)
  expect_invisible(EpiStrainDynamics:::validate_positive_whole_number(10, "test_param"))
})

test_that("validate_positive_whole_number() rejects non-numeric inputs", {
  expect_error(EpiStrainDynamics:::validate_positive_whole_number("5", "test_param"),
               "Argument test_param must be numeric")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(TRUE, "test_param"),
               "Argument test_param must be numeric")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(list(5), "test_param"),
               "Argument test_param must be numeric")
})

test_that("validate_positive_whole_number() rejects multiple values", {
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(c(1, 2, 3), "test_param"),
               "Argument test_param must be a single value")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(c(5, 10), "test_param"),
               "Argument test_param must be a single value")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(integer(0), "test_param"),
               "Argument test_param must be a single value")
})

test_that("validate_positive_whole_number() rejects non-finite values", {
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(NA, "test_param"),
               "Argument test_param must be numeric")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(Inf, "test_param"),
               "Argument test_param must be a finite number")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(-Inf, "test_param"),
               "Argument test_param must be a finite number")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(NaN, "test_param"),
               "Argument test_param must be a finite number")
})

test_that("validate_positive_whole_number() rejects non-whole numbers", {
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(1.5, "test_param"),
               "Argument test_param must be a whole number")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(3.14, "test_param"),
               "Argument test_param must be a whole number")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(-2.7, "test_param"),
               "Argument test_param must be a whole number")
})

test_that("validate_positive_whole_number() rejects non-positive numbers", {
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(0, "test_param"),
               "Argument test_param must be a positive number")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(-1, "test_param"),
               "Argument test_param must be a positive number")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(-10, "test_param"),
               "Argument test_param must be a positive number")
})

test_that("validate_positive_whole_number() uses correct argument name in errors", {
  expect_error(EpiStrainDynamics:::validate_positive_whole_number(-1, "spline_degree"),
               "Argument spline_degree must be a positive number")
  expect_error(EpiStrainDynamics:::validate_positive_whole_number("abc", "days_per_knot"),
               "Argument days_per_knot must be numeric")
})

# Tests for validate_class_inherits()
test_that("validate_class_inherits() passes valid class", {
  df <- data.frame(x = 1:3, y = letters[1:3])
  mat <- matrix(1:6, nrow = 2)

  expect_silent(EpiStrainDynamics:::validate_class_inherits(df, "data.frame"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(mat, "matrix"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(1:5, "integer"))
})

test_that("validate_class_inherits() handles inheritance", {
  # Create object with multiple classes
  obj <- structure(list(a = 1), class = c("child_class", "parent_class"))

  expect_silent(EpiStrainDynamics:::validate_class_inherits(obj, "parent_class"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(obj, "child_class"))
})

# Tests for validate_gi_dist()
test_that("validate_gi_dist accepts valid functions", {
  # Standard valid function
  gi_dist <- function(x) 4 * x * exp(-2 * x)
  expect_silent(validate_gi_dist(gi_dist))

  # Exponential-like
  gi_dist1 <- function(x) exp(-x)
  expect_silent(validate_gi_dist(gi_dist1))

  # Gamma-like
  gi_dist2 <- function(x) x^2 * exp(-x)
  expect_silent(validate_gi_dist(gi_dist2))

  # Constant
  gi_dist3 <- function(x) rep(0.1, length(x))
  expect_silent(validate_gi_dist(gi_dist3))
})

test_that("validate_gi_dist rejects invalid input types", {
  # Non-function inputs
  expect_error(validate_gi_dist(5), class = "rlang_error")
  expect_error(validate_gi_dist("not a function"), class = "rlang_error")

  # Function with no arguments
  gi_dist <- function() 1
  expect_error(validate_gi_dist(gi_dist), class = "rlang_error")
})

test_that("validate_gi_dist rejects functions with invalid outputs", {
  # Function that errors
  gi_dist1 <- function(x) stop("intentional error")
  expect_error(validate_gi_dist(gi_dist1), class = "rlang_error")

  # Non-numeric output
  gi_dist2 <- function(x) as.character(x)
  expect_error(validate_gi_dist(gi_dist2), class = "rlang_error")

  # Non-vectorized function
  gi_dist3 <- function(x) sum(x)
  expect_error(validate_gi_dist(gi_dist3), class = "rlang_error")

  # Negative values
  gi_dist4 <- function(x) -x
  expect_error(validate_gi_dist(gi_dist4), class = "rlang_error")
})

test_that("validate_gi_dist warns about non-finite values", {
  gi_dist <- function(x) ifelse(x == 0, Inf, x)
  expect_warning(validate_gi_dist(gi_dist), class = "rlang_warning")
})
test_that("validate_class_inherits() rejects wrong class", {
  expect_error(EpiStrainDynamics:::validate_class_inherits("string", "data.frame"),
               "Input must be of class data.frame but got class: character")
  expect_error(EpiStrainDynamics:::validate_class_inherits(1:5, "matrix"),
               "Input must be of class matrix but got class: integer")
  expect_error(EpiStrainDynamics:::validate_class_inherits(list(a = 1), "data.frame"),
               "Input must be of class data.frame but got class: list")
})

test_that("validate_pathogen_combination accepts valid inputs", {
  pathogen_names <- c("alpha", "delta", "omicron", "other")

  # NULL is valid for both
  expect_silent(validate_pathogen_combination(NULL, pathogen_names, "numerator_combination"))
  expect_silent(validate_pathogen_combination(NULL, pathogen_names, "denominator_combination"))

  # Single and multiple pathogens valid for both
  expect_silent(validate_pathogen_combination("alpha", pathogen_names, "numerator_combination"))
  expect_silent(validate_pathogen_combination(c("alpha", "delta"), pathogen_names, "denominator_combination"))

})

test_that("validate_pathogen_combination rejects invalid inputs", {
  pathogen_names <- c("alpha", "delta", "omicron", "other")

  # Non-character input
  expect_error(
    validate_pathogen_combination(123, pathogen_names, "numerator_combination"),
    class = "rlang_error"
  )

  # Invalid pathogen names
  expect_error(
    validate_pathogen_combination("beta", pathogen_names, "denominator_combination"),
    class = "rlang_error"
  )

})
