# Tests for case-insensitive argument matching functions

# Tests for match_arg_case_insensitive()
test_that("match_arg_case_insensitive() handles basic matching correctly", {
  choices <- c("apple", "banana", "cherry")

  # Test exact lowercase match
  expect_equal(match_arg_case_insensitive("apple", choices), "apple")
  expect_equal(match_arg_case_insensitive("banana", choices), "banana")

  # Test uppercase input
  expect_equal(match_arg_case_insensitive("APPLE", choices), "apple")
  expect_equal(match_arg_case_insensitive("BANANA", choices), "banana")

  # Test mixed case input
  expect_equal(match_arg_case_insensitive("Apple", choices), "apple")
  expect_equal(match_arg_case_insensitive("BaNaNa", choices), "banana")
  expect_equal(match_arg_case_insensitive("ChErRy", choices), "cherry")
})

test_that("match_arg_case_insensitive() handles partial matching", {
  choices <- c("independent", "correlated", "shared")

  # Test partial matches that are unambiguous
  expect_equal(match_arg_case_insensitive("ind", choices), "independent")
  expect_equal(match_arg_case_insensitive("cor", choices), "correlated")
  expect_equal(match_arg_case_insensitive("sha", choices), "shared")

  # Test partial matches with different cases
  expect_equal(match_arg_case_insensitive("IND", choices), "independent")
  expect_equal(match_arg_case_insensitive("Cor", choices), "correlated")
})

test_that("match_arg_case_insensitive() handles missing argument", {
  choices <- c("first", "second", "third")

  # Missing argument should return first choice
  result <- match_arg_case_insensitive(choices = choices)
  expect_equal(result, "first")
})

test_that("match_arg_case_insensitive() handles default parameter pattern", {
  choices <- c("apple", "banana", "cherry")

  # Simulate the pattern where arg is the full choices vector (R's default behavior)
  result <- match_arg_case_insensitive(choices, choices)
  expect_equal(result, "apple")

  # Test with mixed case in the full vector
  mixed_choices <- c("Apple", "BANANA", "cherry")
  result <- match_arg_case_insensitive(mixed_choices, choices)
  expect_equal(result, "apple")
})

test_that("match_arg_case_insensitive() handles multiple matches with several.ok", {
  choices <- c("apple", "apricot", "banana")

  # Test multiple matches not allowed by default
  expect_error(
    match_arg_case_insensitive(c("apple", "banana"), choices),
    "argument must be of length 1"
  )

  # Test multiple matches allowed when several.ok = TRUE
  result <- match_arg_case_insensitive(c("apple", "banana"), choices, several.ok = TRUE)
  expect_equal(result, c("apple", "banana"))

  # Test with mixed case
  result <- match_arg_case_insensitive(c("APPLE", "Banana"), choices, several.ok = TRUE)
  expect_equal(result, c("apple", "banana"))
})

test_that("match_arg_case_insensitive() provides informative error messages", {
  choices <- c("red", "green", "blue")

  # Test invalid argument without arg_name
  expect_error(
    match_arg_case_insensitive("yellow", choices),
    "argument should be one of: `red`, `green`, and `blue`. Got `yellow`."
  )

  # Test invalid argument with arg_name
  expect_error(
    match_arg_case_insensitive("yellow", choices, arg_name = "color"),
    "'color' should be one of: `red`, `green`, and `blue`. Got `yellow`."
  )

  # Test multiple invalid arguments
  expect_error(
    match_arg_case_insensitive(c("yellow", "purple"), choices, several.ok = TRUE),
    "argument should be one of: `red`, `green`, and `blue`. Got `yellow` and `purple`."
  )

  # Test length error with arg_name
  expect_error(
    match_arg_case_insensitive(c("red", "green"), choices, arg_name = "color"),
    "'color' must be of length 1"
  )
})

test_that("match_arg_case_insensitive() handles edge cases", {
  choices <- c("a", "b", "c")

  # Test empty string
  expect_error(
    match_arg_case_insensitive("", choices),
    "argument should be one of"
  )

  # Test single character choices
  expect_equal(match_arg_case_insensitive("A", choices), "a")
  expect_equal(match_arg_case_insensitive("b", choices), "b")

  # Test choices with special characters
  special_choices <- c("test-1", "test_2", "test.3")
  expect_equal(match_arg_case_insensitive("TEST-1", special_choices), "test-1")
  expect_equal(match_arg_case_insensitive("test_2", special_choices), "test_2")
})

test_that("match_arg_case_insensitive() preserves original choice case in output", {
  choices <- c("CamelCase", "snake_case", "UPPER_CASE")

  # Input case should not affect output case - should return original choice case
  expect_equal(match_arg_case_insensitive("camelcase", choices), "CamelCase")
  expect_equal(match_arg_case_insensitive("SNAKE_CASE", choices), "snake_case")
  expect_equal(match_arg_case_insensitive("upper_case", choices), "UPPER_CASE")
})

# Tests for match_args_case_insensitive()
test_that("match_args_case_insensitive() handles basic multiple argument matching", {
  arg_list <- list(
    color = "RED",
    size = "large",
    shape = "Circle"
  )

  choices_list <- list(
    color = c("red", "green", "blue"),
    size = c("small", "medium", "large"),
    shape = c("circle", "square", "triangle")
  )

  result <- match_args_case_insensitive(arg_list, choices_list)

  expected <- list(
    color = "red",
    size = "large",
    shape = "circle"
  )

  expect_equal(result, expected)
})

test_that("match_args_case_insensitive() handles partial matching for multiple args", {
  arg_list <- list(
    method = "rand",
    structure = "indep"
  )

  choices_list <- list(
    method = c("random-walk", "p-spline"),
    structure = c("independent", "correlated", "shared")
  )

  result <- match_args_case_insensitive(arg_list, choices_list)

  expected <- list(
    method = "random-walk",
    structure = "independent"
  )

  expect_equal(result, expected)
})

test_that("match_args_case_insensitive() handles several.ok parameter", {
  arg_list <- list(
    colors = c("red", "blue"),
    sizes = c("small", "large")
  )

  choices_list <- list(
    colors = c("red", "green", "blue"),
    sizes = c("small", "medium", "large")
  )

  # Should fail with several.ok = FALSE (default)
  expect_error(
    match_args_case_insensitive(arg_list, choices_list),
    "'colors' must be of length 1"
  )

  # Should work with several.ok = TRUE
  result <- match_args_case_insensitive(arg_list, choices_list, several.ok = TRUE)
  expected <- list(
    colors = c("red", "blue"),
    sizes = c("small", "large")
  )

  expect_equal(result, expected)
})

test_that("match_args_case_insensitive() validates input arguments", {
  # Test non-list inputs
  expect_error(
    match_args_case_insensitive("not_a_list", list(a = c("x", "y"))),
    "Both `arg_list` and `choices_list` must be lists"
  )

  expect_error(
    match_args_case_insensitive(list(a = "x"), "not_a_list"),
    "Both `arg_list` and `choices_list` must be lists"
  )

  # Test missing choices for arguments
  arg_list <- list(color = "red", size = "large")
  choices_list <- list(color = c("red", "blue"))  # missing 'size'

  expect_error(
    match_args_case_insensitive(arg_list, choices_list),
    "Missing choices for arguments: size"
  )

  # Test multiple missing choices
  arg_list <- list(color = "red", size = "large", shape = "circle")
  choices_list <- list(color = c("red", "blue"))  # missing 'size' and 'shape'

  expect_error(
    match_args_case_insensitive(arg_list, choices_list),
    "Missing choices for arguments: size and shape"
  )
})

test_that("match_args_case_insensitive() propagates individual matching errors", {
  arg_list <- list(
    color = "yellow",  # invalid choice
    size = "large"     # valid choice
  )

  choices_list <- list(
    color = c("red", "green", "blue"),
    size = c("small", "medium", "large")
  )

  # Should get error from the invalid color choice
  expect_error(
    match_args_case_insensitive(arg_list, choices_list),
    "'color' should be one of: `red`, `green`, and `blue`. Got `yellow`."
  )
})

test_that("match_args_case_insensitive() handles empty lists", {
  # Empty arg_list should return empty result
  result <- match_args_case_insensitive(list(), list())
  expect_equal(result, list())
  expect_length(result, 0)

  # Non-empty choices_list with empty arg_list should work
  result <- match_args_case_insensitive(list(), list(color = c("red", "blue")))
  expect_equal(result, list())
})

test_that("match_args_case_insensitive() preserves argument order", {
  arg_list <- list(
    z_param = "option1",
    a_param = "option2",
    m_param = "option3"
  )

  choices_list <- list(
    z_param = c("option1", "option2"),
    a_param = c("option1", "option2"),
    m_param = c("option1", "option2", "option3")
  )

  result <- match_args_case_insensitive(arg_list, choices_list)

  # Names should preserve original order
  expect_equal(names(result), c("z_param", "a_param", "m_param"))
  expect_equal(result$z_param, "option1")
  expect_equal(result$a_param, "option2")
  expect_equal(result$m_param, "option3")
})

# Integration tests using both functions together
test_that("functions work together for common use cases", {
  # Simulate a function that uses both single and multiple argument matching

  # Single argument case
  smoothing_choice <- match_arg_case_insensitive(
    "INDEPENDENT",
    c("shared", "independent", "correlated"),
    arg_name = "smoothing_structure"
  )
  expect_equal(smoothing_choice, "independent")

  # Multiple argument case
  model_params <- list(
    smoothing = "CORR",
    noise = "PATH"
  )

  param_choices <- list(
    smoothing = c("shared", "independent", "correlated"),
    noise = c("observation_noise_only", "pathogen_specific_noise")
  )

  matched_params <- match_args_case_insensitive(model_params, param_choices)

  expect_equal(matched_params$smoothing, "correlated")
  expect_equal(matched_params$noise, "pathogen_specific_noise")
})

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

# Tests for validate_list_vector_lengths()
test_that("validate_list_vector_lengths() passes valid lists", {
  valid_list <- list(a = 1:3, b = letters[1:3], c = c(1.1, 2.2, 3.3))
  expect_silent(validate_list_vector_lengths(valid_list, "test_list"))

  single_element <- list(x = 1:5)
  expect_silent(validate_list_vector_lengths(single_element, "single_list"))
})

test_that("validate_list_vector_lengths() validates against reference length", {
  valid_list <- list(a = 1:3, b = letters[1:3])
  expect_silent(validate_list_vector_lengths(valid_list, "test_list", reference_length = 3))

  expect_error(validate_list_vector_lengths(valid_list, "test_list", reference_length = 5),
               "Vectors in test_list must have length 5 but found length 3")
})

test_that("validate_list_vector_lengths() rejects empty lists", {
  empty_list <- list()
  expect_error(validate_list_vector_lengths(empty_list, "empty_list"),
               "empty_list cannot be empty")
})

test_that("validate_list_vector_lengths() rejects mismatched lengths", {
  invalid_list <- list(short = 1:2, long = 1:5)
  expect_error(validate_list_vector_lengths(invalid_list, "mixed_list"),
               "All vectors in mixed_list must have the same length")
})

test_that("validate_list_vector_lengths() includes detailed error information", {
  invalid_list <- list(vec1 = 1:2, vec2 = 1:4, vec3 = 1:3)
  expect_error(validate_list_vector_lengths(invalid_list, "test_list"),
               "vec1 \\(length 2 \\).*vec2 \\(length 4 \\).*vec3 \\(length 3 \\)")
})

test_that("validate_list_vector_lengths() uses correct list name in errors", {
  invalid_list <- list(a = 1:2, b = 1:3)
  expect_error(validate_list_vector_lengths(invalid_list, "my_data"),
               "All vectors in my_data must have the same length")

  empty_list <- list()
  expect_error(validate_list_vector_lengths(empty_list, "strain_data"),
               "strain_data cannot be empty")
})

