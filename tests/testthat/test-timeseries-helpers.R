# Tests for timeseries helper functions
# These tests focus on validation checks and edge cases in the helper functions

# ==============================================================================
# Tests for is_timeseries_class()
# ==============================================================================

test_that("is_timeseries_class() correctly identifies time series objects", {
  # Test with data.frame (should be FALSE)
  df <- data.frame(x = 1:10, y = 1:10)
  expect_false(is_timeseries_class(df))

  # Test with tibble (should be FALSE)
  if (requireNamespace("tibble", quietly = TRUE)) {
    tbl <- tibble::tibble(x = 1:10, y = 1:10)
    expect_false(is_timeseries_class(tbl))
  }

  # Test with ts (should be TRUE)
  ts_obj <- ts(1:10, start = 1, frequency = 1)
  expect_true(is_timeseries_class(ts_obj))

  # Test with mts (should be TRUE)
  mts_obj <- ts(matrix(1:20, ncol = 2), start = 1, frequency = 1)
  expect_true(is_timeseries_class(mts_obj))

  # Test with xts (should be TRUE)
  if (requireNamespace("xts", quietly = TRUE)) {
    xts_obj <- xts::xts(1:10, order.by = Sys.Date() + 0:9)
    expect_true(is_timeseries_class(xts_obj))
  }

  # Test with zoo (should be TRUE)
  if (requireNamespace("zoo", quietly = TRUE)) {
    zoo_obj <- zoo::zoo(1:10, order.by = 1:10)
    expect_true(is_timeseries_class(zoo_obj))
  }

  # Test with tsibble (should be TRUE)
  if (requireNamespace("tsibble", quietly = TRUE)) {
    tsbl_obj <- tsibble::tsibble(
      time = 1:10,
      value = 1:10,
      index = time
    )
    expect_true(is_timeseries_class(tsbl_obj))
  }
})

# ==============================================================================
# Tests for convert_ts_to_tsibble()
# ==============================================================================

test_that("convert_ts_to_tsibble() handles xts objects correctly", {
  skip_if_not_installed("xts")
  skip_if_not_installed("timetk")

  # Create xts object
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10)
  xts_obj <- xts::xts(
    data.frame(a = 1:10, b = 11:20),
    order.by = dates
  )

  result <- convert_ts_to_tsibble(xts_obj)

  expect_s3_class(result, "tbl_ts")
  expect_true("index" %in% names(result))
  expect_equal(nrow(result), 10)
})

test_that("convert_ts_to_tsibble() handles zoo objects correctly", {
  skip_if_not_installed("zoo")
  skip_if_not_installed("timetk")

  # Create zoo object
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10)
  zoo_obj <- zoo::zoo(
    data.frame(a = 1:10, b = 11:20),
    order.by = dates
  )

  result <- convert_ts_to_tsibble(zoo_obj)

  expect_s3_class(result, "tbl_ts")
  expect_true("index" %in% names(result))
  expect_equal(nrow(result), 10)
})

test_that("convert_ts_to_tsibble() keeps multivariate ts in wide format", {
  # Create multivariate time series
  mts_obj <- ts(matrix(1:30, ncol = 3), start = 1, frequency = 1)
  colnames(mts_obj) <- c("a", "b", "c")

  result <- convert_ts_to_tsibble(mts_obj)

  expect_s3_class(result, "tbl_ts")
  # Should have 3 value columns plus 1 index column
  expect_true(all(c("a", "b", "c") %in% names(result)))
  expect_equal(nrow(result), 10)
})

test_that("convert_ts_to_tsibble() handles univariate ts correctly", {
  # Create univariate time series
  ts_obj <- ts(1:10, start = 1, frequency = 1)

  result <- convert_ts_to_tsibble(ts_obj)

  expect_s3_class(result, "tbl_ts")
  expect_equal(nrow(result), 10)
})

test_that("convert_ts_to_tsibble() passes through tsibble objects", {
  skip_if_not_installed("tsibble")

  # Create tsibble
  tsbl_obj <- tsibble::tsibble(
    time = 1:10,
    value = 1:10,
    index = time
  )

  result <- convert_ts_to_tsibble(tsbl_obj)

  expect_identical(result, tsbl_obj)
})

test_that("convert_ts_to_tsibble() errors on unsupported classes", {
  # Create an unsupported object
  unsupported <- list(x = 1:10)
  class(unsupported) <- "unsupported_ts_class"

  expect_error(
    convert_ts_to_tsibble(unsupported),
    "Unable to convert time series object"
  )
})

# ==============================================================================
# Tests for create_validated_timeseries() - tsibble validation checks
# ==============================================================================

test_that("create_validated_timeseries() detects time series with gaps", {
  # Create data with gaps
  df <- data.frame(
    date = as.Date(c("2020-01-01", "2020-01-02", "2020-01-04", "2020-01-05")), # Missing day 3
    cases = c(10, 15, 20, 25),
    alpha = c(5, 8, 10, 12)
  )

  expect_error(
    create_validated_timeseries(
      data = df,
      columns = c("cases", "alpha"),
      time_col = "date"
    ),
    "Time series has gaps"
  )
})

test_that("create_validated_timeseries() detects irregular time series", {
  # Create irregularly spaced data
  df <- data.frame(
    date = as.Date(c("2020-01-01", "2020-01-02", "2020-01-04", "2020-01-08")), # Irregular spacing
    cases = c(10, 15, 20, 25),
    alpha = c(5, 8, 10, 12)
  )

  expect_error(
    create_validated_timeseries(
      data = df,
      columns = c("cases", "alpha"),
      time_col = "date"
    ),
    "Time series has gaps"
  )
})

test_that("create_validated_timeseries() detects unordered time series", {
  # Create unordered data
  df <- data.frame(
    date = as.Date(c("2020-01-01", "2020-01-03", "2020-01-02", "2020-01-04")), # Out of order
    cases = c(10, 15, 20, 25),
    alpha = c(5, 8, 10, 12)
  )

  expect_error(
    create_validated_timeseries(
      data = df,
      columns = c("cases", "alpha"),
      time_col = "date"
    ),
    "Time series is not ordered"
  )
})

test_that("create_validated_timeseries() warns when time arg provided for ts objects", {
  skip_if_not_installed("xts")
  skip_if_not_installed("timetk")

  # Create xts object
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10)
  xts_obj <- xts::xts(
    data.frame(cases = 1:10, alpha = 11:20),
    order.by = dates
  )

  expect_warning(
    create_validated_timeseries(
      data = xts_obj,
      columns = c("cases", "alpha"),
      time_col = "date"  # This should trigger a warning
    ),
    "time.*argument is ignored when data is a time series object"
  )
})

test_that("create_validated_timeseries() returns valid tsibble for correct data", {
  # Create valid data
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = 1:10,
    alpha = 11:20
  )

  result <- create_validated_timeseries(
    data = df,
    columns = c("cases", "alpha"),
    time_col = "date"
  )

  expect_s3_class(result, "tbl_ts")
  expect_equal(nrow(result), 10)
  expect_true(all(c("date", "cases", "alpha") %in% names(result)))
  expect_false(tsibble::has_gaps(result)$.gaps)
  expect_true(tsibble::is_regular(result))
  expect_true(tsibble::is_ordered(result))
})

test_that("create_validated_timeseries() works with time series objects", {
  skip_if_not_installed("xts")
  skip_if_not_installed("timetk")

  # Create valid xts object
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10)
  xts_obj <- xts::xts(
    data.frame(cases = 1:10, alpha = 11:20),
    order.by = dates
  )

  result <- create_validated_timeseries(
    data = xts_obj,
    columns = c("cases", "alpha"),
    time_col = NULL
  )

  expect_s3_class(result, "tbl_ts")
  expect_equal(nrow(result), 10)
  expect_true(all(c("cases", "alpha") %in% names(result)))
  expect_false(tsibble::has_gaps(result)$.gaps)
  expect_true(tsibble::is_regular(result))
  expect_true(tsibble::is_ordered(result))
})

test_that("create_validated_timeseries() handles data.table correctly", {
  skip_if_not_installed("data.table")

  # Create valid data.table
  dt <- data.table::data.table(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = 1:10,
    alpha = 11:20
  )

  result <- create_validated_timeseries(
    data = dt,
    columns = c("cases", "alpha"),
    time_col = "date"
  )

  expect_s3_class(result, "tbl_ts")
  expect_equal(nrow(result), 10)
  expect_true(all(c("date", "cases", "alpha") %in% names(result)))
})

test_that("create_validated_timeseries() handles tibble correctly", {
  skip_if_not_installed("tibble")

  # Create valid tibble
  tbl <- tibble::tibble(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = 1:10,
    alpha = 11:20
  )

  result <- create_validated_timeseries(
    data = tbl,
    columns = c("cases", "alpha"),
    time_col = "date"
  )

  expect_s3_class(result, "tbl_ts")
  expect_equal(nrow(result), 10)
  expect_true(all(c("date", "cases", "alpha") %in% names(result)))
})

test_that("create_validated_timeseries() subsets to only required columns", {
  # Create data with extra columns
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = 1:10,
    alpha = 11:20,
    extra1 = 21:30,
    extra2 = 31:40
  )

  result <- create_validated_timeseries(
    data = df,
    columns = c("cases", "alpha"),
    time_col = "date"
  )

  expect_s3_class(result, "tbl_ts")
  expect_equal(ncol(result), 3)  # date + cases + alpha
  expect_true(all(c("date", "cases", "alpha") %in% names(result)))
  expect_false("extra1" %in% names(result))
  expect_false("extra2" %in% names(result))
})

# ==============================================================================
# Tests for error handling in tsibble creation
# ==============================================================================

test_that("create_validated_timeseries() handles duplicate time indices", {
  # Create data with duplicate dates
  df <- data.frame(
    date = as.Date(c("2020-01-01", "2020-01-02", "2020-01-02", "2020-01-03")), # Duplicate
    cases = 1:4,
    alpha = 5:8
  )

  expect_error(
    create_validated_timeseries(
      data = df,
      columns = c("cases", "alpha"),
      time_col = "date"
    ),
    "Error creating tsibble|duplicated"
  )
})

# ==============================================================================
# Tests for numeric column validation
# ==============================================================================

test_that("create_validated_timeseries() validates numeric columns for data frames", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = 1:10,
    alpha = as.character(11:20)  # Non-numeric
  )

  expect_error(
    create_validated_timeseries(
      data = df,
      columns = c("cases", "alpha"),
      time_col = "date"
    ),
    "numeric"
  )
})

test_that("create_validated_timeseries() validates numeric columns for tsibble objects", {
  skip_if_not_installed("tsibble")

  # Create tsibble with character column
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10)
  tsbl_obj <- tsibble::tsibble(
    date = dates,
    cases = 1:10,
    alpha = as.character(11:20),  # Non-numeric
    index = date
  )

  expect_error(
    create_validated_timeseries(
      data = tsbl_obj,
      columns = c("cases", "alpha"),
      time_col = NULL
    ),
    "numeric"
  )
})

# ==============================================================================
# Integration tests with actual intake functions
# ==============================================================================

test_that("intake functions fail with gaps in time series", {
  check_package_data()

  # Create data with gaps
  test_data <- sarscov2[c(1:10, 12:20), ]  # Skip row 11

  expect_error(
    single(
      data = test_data,
      case_timeseries = "cases",
      time = "date"
    ),
    "Time series has gaps"
  )
})

test_that("intake functions fail with unordered time series", {
  check_package_data()

  # Create unordered data
  test_data <- sarscov2[c(1:5, 10:6, 11:15), ]  # Swap order

  expect_error(
    single(
      data = test_data,
      case_timeseries = "cases",
      time = "date"
    ),
    "Time series is not ordered"
  )
})

test_that("multiple() fails with duplicate time indices", {
  check_package_data()

  # Create data with duplicate dates
  test_data <- rbind(sarscov2[1:10, ], sarscov2[10, ])  # Duplicate last row

  expect_error(
    multiple(
      data = test_data,
      case_timeseries = "cases",
      time = "date",
      component_pathogen_timeseries = c("alpha", "delta")
    ),
    "Error creating tsibble|duplicated"
  )
})

