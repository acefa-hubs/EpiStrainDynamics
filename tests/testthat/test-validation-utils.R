#' @srrstats {G5.2, G5.2b} error and warning behaviour explicitly tested,
#'   including conditions that trigger these messages and comparing with
#'   expected results
#' @srrstats {G3.0} comparisons made between appropriate values

# ==============================================================================
# TESTS: validate_positive_whole_number()
# ==============================================================================

test_that("validate_positive_whole_number() passes valid inputs", {
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(1, "test_param"))
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(5L, "test_param"))
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(100, "test_param"))
  expect_silent(EpiStrainDynamics:::validate_positive_whole_number(3.0, "test_param"))
})

test_that("validate_positive_whole_number() returns invisible NULL", {
  result <- EpiStrainDynamics:::validate_positive_whole_number(5, "test_param")
  expect_null(result)
  expect_invisible(EpiStrainDynamics:::validate_positive_whole_number(10, "test_param"))
})

test_that("validate_positive_whole_number() rejects non-numeric inputs", {
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number("5", "test_param"),
    "Argument test_param must be numeric"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(TRUE, "test_param"),
    "Argument test_param must be numeric"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(list(5), "test_param"),
    "Argument test_param must be numeric"
  )
})

test_that("validate_positive_whole_number() rejects multiple values", {
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(c(1, 2, 3), "test_param"),
    "Argument test_param must be a single value"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(c(5, 10), "test_param"),
    "Argument test_param must be a single value"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(integer(0), "test_param"),
    "Argument test_param must be a single value"
  )
})

test_that("validate_positive_whole_number() rejects non-finite values", {
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(NA, "test_param"),
    "Argument test_param must be numeric"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(Inf, "test_param"),
    "Argument test_param must be a finite number"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(-Inf, "test_param"),
    "Argument test_param must be a finite number"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(NaN, "test_param"),
    "Argument test_param must be a finite number"
  )
})

test_that("validate_positive_whole_number() rejects non-whole numbers", {
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(1.5, "test_param"),
    "Argument test_param must be a whole number"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(3.14, "test_param"),
    "Argument test_param must be a whole number"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(-2.7, "test_param"),
    "Argument test_param must be a whole number"
  )
})

test_that("validate_positive_whole_number() rejects non-positive numbers", {
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(0, "test_param"),
    "Argument test_param must be a positive number"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(-1, "test_param"),
    "Argument test_param must be a positive number"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(-10, "test_param"),
    "Argument test_param must be a positive number"
  )
})

test_that("validate_positive_whole_number() uses correct argument name in errors", {
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number(-1, "spline_degree"),
    "Argument spline_degree must be a positive number"
  )
  expect_error(
    EpiStrainDynamics:::validate_positive_whole_number("abc", "days_per_knot"),
    "Argument days_per_knot must be numeric"
  )
})

# ==============================================================================
# TESTS: validate_class_inherits()
# ==============================================================================

test_that("validate_class_inherits() passes valid class", {
  df <- data.frame(x = 1:3, y = letters[1:3])
  mat <- matrix(1:6, nrow = 2)

  expect_silent(EpiStrainDynamics:::validate_class_inherits(df, "data.frame"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(mat, "matrix"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(1:5, "integer"))
})

test_that("validate_class_inherits() handles inheritance", {
  obj <- structure(list(a = 1), class = c("child_class", "parent_class"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(obj, "parent_class"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(obj, "child_class"))
})

test_that("validate_class_inherits() rejects wrong class", {
  expect_error(
    EpiStrainDynamics:::validate_class_inherits("string", "data.frame"),
    "Input must be of class data.frame but got class: character"
  )
  expect_error(
    EpiStrainDynamics:::validate_class_inherits(1:5, "matrix"),
    "Input must be of class matrix but got class: integer"
  )
  expect_error(
    EpiStrainDynamics:::validate_class_inherits(list(a = 1), "data.frame"),
    "Input must be of class data.frame but got class: list"
  )
})

test_that("validate_class_inherits() handles require_all = FALSE", {
  obj <- structure(list(a = 1), class = c("classA", "classB"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(obj,
    c("classA", "classC"),
    require_all = FALSE
  ))
  expect_error(
    EpiStrainDynamics:::validate_class_inherits(obj,
      c("classX", "classY"),
      require_all = FALSE
    ),
    "Input must inherit from at least one of"
  )
})

test_that("validate_class_inherits() handles multiple required classes", {
  obj <- structure(list(a = 1), class = c("classA", "classB"))
  expect_silent(EpiStrainDynamics:::validate_class_inherits(obj, c("classA", "classB")))
  expect_error(
    EpiStrainDynamics:::validate_class_inherits(obj, c("classA", "classC")),
    "Input must inherit from all classes"
  )
})

test_that("validate_class_inherits() errors on empty class_names", {
  expect_error(
    EpiStrainDynamics:::validate_class_inherits(1:5, character(0)),
    "must be a non-empty character vector"
  )
})

# ==============================================================================
# TESTS: validate_priors()
# ==============================================================================

test_that("validate_priors() accepts valid priors", {
  result <- EpiStrainDynamics:::validate_priors(mean = 0.5, sd = 1.0)
  expect_s3_class(result, "EpiStrainDynamics.prior")
  expect_equal(result$mean, 0.5)
  expect_equal(result$sd, 1.0)
})

test_that("validate_priors() accepts vector priors", {
  result <- EpiStrainDynamics:::validate_priors(mean = c(0.1, 0.2), sd = c(1.0, 1.0))
  expect_s3_class(result, "EpiStrainDynamics.prior")
  expect_length(result$mean, 2)
})

test_that("validate_priors() accepts NULL priors", {
  result <- EpiStrainDynamics:::validate_priors(mean = NULL, sd = NULL)
  expect_s3_class(result, "EpiStrainDynamics.prior")
  expect_null(result$mean)
  expect_null(result$sd)
})

test_that("validate_priors() errors when only one of mean/sd provided", {
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = 0.5, sd = NULL),
    "both mean and sd must be provided"
  )
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = NULL, sd = 1.0),
    "both mean and sd must be provided"
  )
})

test_that("validate_priors() errors on non-numeric priors", {
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = "a", sd = 1.0),
    "must be numeric"
  )
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = 0.5, sd = "b"),
    "must be numeric"
  )
})

test_that("validate_priors() errors on NA values", {
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = NA_real_, sd = 1.0),
    "cannot contain NA values"
  )
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = 0.5, sd = NA_real_),
    "cannot contain NA values"
  )
})

test_that("validate_priors() errors when mean and sd have different lengths", {
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = c(0.1, 0.2), sd = 1.0),
    "must have the same length"
  )
})

test_that("validate_priors() errors on negative values", {
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = -0.5, sd = 1.0),
    "must be positive"
  )
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = 0.5, sd = -1.0),
    "must be positive"
  )
  expect_error(
    EpiStrainDynamics:::validate_priors(mean = 0.5, sd = 0),
    "must be positive"
  )
})

# ==============================================================================
# TESTS: validate_gi_dist()
# ==============================================================================

test_that("validate_gi_dist accepts valid functions", {
  gi_dist <- function(x) 4 * x * exp(-2 * x)
  expect_silent(EpiStrainDynamics:::validate_gi_dist(gi_dist))

  gi_dist1 <- function(x) exp(-x)
  expect_silent(EpiStrainDynamics:::validate_gi_dist(gi_dist1))

  gi_dist2 <- function(x) x^2 * exp(-x)
  expect_silent(EpiStrainDynamics:::validate_gi_dist(gi_dist2))

  gi_dist3 <- function(x) rep(0.1, length(x))
  expect_silent(EpiStrainDynamics:::validate_gi_dist(gi_dist3))
})

test_that("validate_gi_dist rejects invalid input types", {
  expect_error(EpiStrainDynamics:::validate_gi_dist(5), class = "rlang_error")
  expect_error(EpiStrainDynamics:::validate_gi_dist("not a function"), class = "rlang_error")

  gi_dist <- function() 1
  expect_error(EpiStrainDynamics:::validate_gi_dist(gi_dist), class = "rlang_error")
})

test_that("validate_gi_dist rejects functions with invalid outputs", {
  gi_dist1 <- function(x) stop("intentional error")
  expect_error(EpiStrainDynamics:::validate_gi_dist(gi_dist1), class = "rlang_error")

  gi_dist2 <- function(x) as.character(x)
  expect_error(EpiStrainDynamics:::validate_gi_dist(gi_dist2), class = "rlang_error")

  gi_dist3 <- function(x) sum(x)
  expect_error(EpiStrainDynamics:::validate_gi_dist(gi_dist3), class = "rlang_error")

  gi_dist4 <- function(x) -x
  expect_error(EpiStrainDynamics:::validate_gi_dist(gi_dist4), class = "rlang_error")
})

test_that("validate_gi_dist warns about non-finite values", {
  gi_dist <- function(x) ifelse(x == 0, Inf, x)
  expect_warning(EpiStrainDynamics:::validate_gi_dist(gi_dist), class = "rlang_warning")
})

# ==============================================================================
# TESTS: validate_pathogen_combination()
# ==============================================================================

test_that("validate_pathogen_combination accepts valid inputs", {
  pathogen_names <- c("alpha", "delta", "omicron", "other")

  expect_silent(EpiStrainDynamics:::validate_pathogen_combination(
    NULL, pathogen_names, "numerator_combination"
  ))
  expect_silent(EpiStrainDynamics:::validate_pathogen_combination(
    NULL, pathogen_names, "denominator_combination"
  ))
  expect_silent(EpiStrainDynamics:::validate_pathogen_combination(
    "alpha", pathogen_names, "numerator_combination"
  ))
  expect_silent(EpiStrainDynamics:::validate_pathogen_combination(
    c("alpha", "delta"), pathogen_names, "denominator_combination"
  ))
})

test_that("validate_pathogen_combination returns correct indices", {
  pathogen_names <- c("alpha", "delta", "omicron", "other")

  result <- EpiStrainDynamics:::validate_pathogen_combination(
    NULL, pathogen_names, "numerator_combination"
  )
  expect_equal(result, 1:4)

  result <- EpiStrainDynamics:::validate_pathogen_combination(
    "delta", pathogen_names, "numerator_combination"
  )
  expect_equal(result, 2L)

  result <- EpiStrainDynamics:::validate_pathogen_combination(
    c("omicron", "alpha"), pathogen_names, "numerator_combination"
  )
  expect_equal(result, c(3L, 1L))
})

test_that("validate_pathogen_combination rejects invalid pathogen names", {
  pathogen_names <- c("alpha", "delta", "omicron", "other")

  expect_error(
    EpiStrainDynamics:::validate_pathogen_combination(
      "beta", pathogen_names, "denominator_combination"
    ),
    class = "rlang_error"
  )
  expect_error(
    EpiStrainDynamics:::validate_pathogen_combination(
      c("alpha", "beta"), pathogen_names, "numerator_combination"
    ),
    class = "rlang_error"
  )
})

# ==============================================================================
# TESTS: check_missing_data()
# ==============================================================================

test_that("check_missing_data() passes when no NAs present", {
  df <- data.frame(cases = c(10, 20, 30), dates = as.Date("2024-01-01") + 0:2)
  expect_silent(EpiStrainDynamics:::check_missing_data(df, "cases", "case_timeseries"))
})

test_that("check_missing_data() errors when NAs present", {
  df <- data.frame(cases = c(10, NA, 30))
  expect_error(
    EpiStrainDynamics:::check_missing_data(df, "cases", "case_timeseries"),
    "Missing values"
  )
})

test_that("check_missing_data() shows correct row indices in error", {
  df <- data.frame(cases = c(10, NA, 30, NA, 50))
  expect_error(
    EpiStrainDynamics:::check_missing_data(df, "cases", "case_timeseries"),
    "row\\(s\\) 2, 4"
  )
})

test_that("check_missing_data() summarises many NAs", {
  df <- data.frame(cases = c(NA, NA, NA, NA, NA, NA, 10))
  expect_error(
    EpiStrainDynamics:::check_missing_data(df, "cases", "case_timeseries"),
    "6 total"
  )
})

# ==============================================================================
# TESTS: check_list_columns()
# ==============================================================================

test_that("check_list_columns() passes for normal data frames", {
  df <- data.frame(cases = 1:5, dates = as.Date("2024-01-01") + 0:4)
  expect_silent(EpiStrainDynamics:::check_list_columns(df, "cases"))
})

test_that("check_list_columns() errors when list columns present", {
  df <- data.frame(x = 1:3)
  df$list_col <- list(1, 2, 3)
  expect_error(
    EpiStrainDynamics:::check_list_columns(df, "list_col"),
    "List columns detected"
  )
})

test_that("check_list_columns() passes when relevant cols don't exist", {
  df <- data.frame(x = 1:3)
  expect_silent(EpiStrainDynamics:::check_list_columns(df, "nonexistent_col"))
})

# ==============================================================================
# TESTS: is_timeseries_class()
# ==============================================================================

test_that("is_timeseries_class() correctly identifies time series objects", {
  expect_true(EpiStrainDynamics:::is_timeseries_class(ts(1:10)))
  expect_false(EpiStrainDynamics:::is_timeseries_class(data.frame(x = 1:5)))
  expect_false(EpiStrainDynamics:::is_timeseries_class(1:10))
  expect_false(EpiStrainDynamics:::is_timeseries_class(list(x = 1:5)))
})
