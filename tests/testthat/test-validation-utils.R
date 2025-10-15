#' @srrstats {G5.2, G5.2b} error and warning behaviour explicitly tested,
#'   including conditions that trigger these messages and comparing with
#'   expected results
#' @srrstats {G3.0} comparisons made between appropriate values
#'
# Tests for validate_positive_whole_number()
test_that("validate_positive_whole_number() passes valid inputs", {
  # Should not throw errors for valid positive whole numbers
  expect_silent(validate_positive_whole_number(1, "test_param"))
  expect_silent(validate_positive_whole_number(5L, "test_param"))
  expect_silent(validate_positive_whole_number(100, "test_param"))
  expect_silent(validate_positive_whole_number(3.0, "test_param"))  # Whole number as double
})

test_that("validate_positive_whole_number() returns invisible NULL", {
  result <- validate_positive_whole_number(5, "test_param")
  expect_null(result)
  expect_invisible(validate_positive_whole_number(10, "test_param"))
})

test_that("validate_positive_whole_number() rejects non-numeric inputs", {
  expect_error(validate_positive_whole_number("5", "test_param"),
               "Argument test_param must be numeric")
  expect_error(validate_positive_whole_number(TRUE, "test_param"),
               "Argument test_param must be numeric")
  expect_error(validate_positive_whole_number(list(5), "test_param"),
               "Argument test_param must be numeric")
})

test_that("validate_positive_whole_number() rejects multiple values", {
  expect_error(validate_positive_whole_number(c(1, 2, 3), "test_param"),
               "Argument test_param must be a single value")
  expect_error(validate_positive_whole_number(c(5, 10), "test_param"),
               "Argument test_param must be a single value")
  expect_error(validate_positive_whole_number(integer(0), "test_param"),
               "Argument test_param must be a single value")
})

test_that("validate_positive_whole_number() rejects non-finite values", {
  expect_error(validate_positive_whole_number(NA, "test_param"),
               "Argument test_param must be numeric")
  expect_error(validate_positive_whole_number(Inf, "test_param"),
               "Argument test_param must be a finite number")
  expect_error(validate_positive_whole_number(-Inf, "test_param"),
               "Argument test_param must be a finite number")
  expect_error(validate_positive_whole_number(NaN, "test_param"),
               "Argument test_param must be a finite number")
})

test_that("validate_positive_whole_number() rejects non-whole numbers", {
  expect_error(validate_positive_whole_number(1.5, "test_param"),
               "Argument test_param must be a whole number")
  expect_error(validate_positive_whole_number(3.14, "test_param"),
               "Argument test_param must be a whole number")
  expect_error(validate_positive_whole_number(-2.7, "test_param"),
               "Argument test_param must be a whole number")
})

test_that("validate_positive_whole_number() rejects non-positive numbers", {
  expect_error(validate_positive_whole_number(0, "test_param"),
               "Argument test_param must be a positive number")
  expect_error(validate_positive_whole_number(-1, "test_param"),
               "Argument test_param must be a positive number")
  expect_error(validate_positive_whole_number(-10, "test_param"),
               "Argument test_param must be a positive number")
})

test_that("validate_positive_whole_number() uses correct argument name in errors", {
  expect_error(validate_positive_whole_number(-1, "spline_degree"),
               "Argument spline_degree must be a positive number")
  expect_error(validate_positive_whole_number("abc", "days_per_knot"),
               "Argument days_per_knot must be numeric")
})

# Tests for validate_class_inherits()
test_that("validate_class_inherits() passes valid class", {
  df <- data.frame(x = 1:3, y = letters[1:3])
  mat <- matrix(1:6, nrow = 2)

  expect_silent(validate_class_inherits(df, "data.frame"))
  expect_silent(validate_class_inherits(mat, "matrix"))
  expect_silent(validate_class_inherits(1:5, "integer"))
})

test_that("validate_class_inherits() handles inheritance", {
  # Create object with multiple classes
  obj <- structure(list(a = 1), class = c("child_class", "parent_class"))

  expect_silent(validate_class_inherits(obj, "parent_class"))
  expect_silent(validate_class_inherits(obj, "child_class"))
})

test_that("validate_class_inherits() rejects wrong class", {
  expect_error(validate_class_inherits("string", "data.frame"),
               "Input must be of class data.frame but got class: character")
  expect_error(validate_class_inherits(1:5, "matrix"),
               "Input must be of class matrix but got class: integer")
  expect_error(validate_class_inherits(list(a = 1), "data.frame"),
               "Input must be of class data.frame but got class: list")
})

# Tests for validate_matching_lengths()
test_that("validate_matching_lengths() passes vectors with same length", {
  expect_silent(validate_matching_lengths(a = 1:3, b = letters[1:3], c = c(1.1, 2.2, 3.3)))
  expect_silent(validate_matching_lengths(x = 1:5, y = 6:10))
  expect_silent(validate_matching_lengths(single = 1))
})

test_that("validate_matching_lengths() handles empty input", {
  expect_silent(validate_matching_lengths())
})

test_that("validate_matching_lengths() rejects different lengths", {
  expect_error(validate_matching_lengths(a = 1:3, b = 1:5),
               "All input vectors must have the same length")
  expect_error(validate_matching_lengths(x = 1:2, y = 1:3, z = 1:4),
               "All input vectors must have the same length")
})

test_that("validate_matching_lengths() includes length information in error", {
  expect_error(validate_matching_lengths(short = 1:2, long = 1:5),
               "short \\(length 2 \\).*long \\(length 5 \\)")
})

# Tests for validate_list_vector()
test_that("validate_list_vector() passes valid lists", {
  valid_list <- list(a = 1:3, b = letters[1:3], c = c(1.1, 2.2, 3.3))
  expect_silent(validate_list_vector(valid_list, "test_list"))

  single_element <- list(x = 1:5)
  expect_silent(validate_list_vector(single_element, "single_list"))
})

test_that("validate_list_vector() validates against reference length", {
  valid_list <- list(a = 1:3, b = letters[1:3])
  expect_silent(validate_list_vector(valid_list, "test_list", reference_length = 3))

  expect_error(validate_list_vector(valid_list, "test_list", reference_length = 5),
               "Vectors in test_list must have length 5 but found length 3")
})

test_that("validate_list_vector() rejects non-lists", {
  expect_error(validate_list_vector("not_a_list", "test_param"),
               "Argument test_param must be a list")
  expect_error(validate_list_vector(c(1, 2, 3), "test_param"),
               "Argument test_param must be a list")
  expect_error(validate_list_vector(data.frame(x = 1:3), "test_param"),
               "Argument test_param must be a list")
})

test_that("validate_list_vector() rejects empty lists", {
  empty_list <- list()
  expect_error(validate_list_vector(empty_list, "empty_list"),
               "empty_list cannot be empty")
})

test_that("validate_list_vector() rejects unnamed lists", {
  unnamed_list <- list(1:3, letters[1:3])
  expect_error(validate_list_vector(unnamed_list, "unnamed_list"),
               "Argument unnamed_list must have named elements")
})

test_that("validate_list_vector() rejects lists with empty or missing names", {
  empty_name_list <- list(1:3, letters[1:3])
  expect_error(validate_list_vector(empty_name_list, "bad_names"),
               "All elements in bad_names must have non-empty names")

  na_name_list <- list(a = 1:3, b = letters[1:3])
  names(na_name_list)[2] <- NA
  expect_error(validate_list_vector(na_name_list, "bad_names"),
               "All elements in bad_names must have non-empty names")
})

test_that("validate_list_vector() rejects mismatched lengths", {
  invalid_list <- list(short = 1:2, long = 1:5)
  expect_error(validate_list_vector(invalid_list, "mixed_list"),
               "All vectors in mixed_list must have the same length")
})

test_that("validate_list_vector() includes detailed error information", {
  invalid_list <- list(vec1 = 1:2, vec2 = 1:4, vec3 = 1:3)
  expect_error(validate_list_vector(invalid_list, "test_list"),
               "vec1 \\(length 2 \\).*vec2 \\(length 4 \\).*vec3 \\(length 3 \\)")
})

test_that("validate_list_vector() uses correct list name in errors", {
  invalid_list <- list(a = 1:2, b = 1:3)
  expect_error(validate_list_vector(invalid_list, "my_data"),
               "All vectors in my_data must have the same length")

  empty_list <- list()
  expect_error(validate_list_vector(empty_list, "strain_data"),
               "strain_data cannot be empty")
})

test_that("validate_list_vector() returns invisible TRUE", {
  valid_list <- list(a = 1:3, b = letters[1:3])
  result <- validate_list_vector(valid_list, "test_list")
  expect_true(result)
  expect_invisible(validate_list_vector(valid_list, "test_list"))
})
