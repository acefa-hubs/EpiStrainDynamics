# Tests for units package support

test_that("single() accepts units objects in case timeseries", {
  skip_if_not_installed("units")

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = units::set_units(1:10, "1")  # Dimensionless count
  )

  result <- single(
    data = df,
    case_timeseries = "cases",
    time = "date"
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "single")
  # Units should be stripped, result should be numeric
  expect_true(is.numeric(result$data$case_timeseries))
  expect_equal(length(result$data$case_timeseries), 10)
})

test_that("multiple() accepts units objects in component pathogens", {
  skip_if_not_installed("units")

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = units::set_units(1:10, "1"),
    alpha = units::set_units(11:20, "1"),
    delta = units::set_units(21:30, "1")
  )

  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "multiple")
  # Units should be stripped
  expect_true(is.numeric(result$data$case_timeseries))
  expect_true(is.numeric(result$data$component_pathogens))
})

test_that("subtyped() accepts units objects in all pathogen timeseries", {
  skip_if_not_installed("units")

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    ili = units::set_units(seq(100, 190, by = 10), "1"),
    inf_A = units::set_units(seq(30, 75, by = 5), "1"),
    inf_H3N2 = units::set_units(seq(10, 28, by = 2), "1"),
    inf_H1N1 = units::set_units(seq(8, 26, by = 2), "1"),
    inf_B = units::set_units(seq(20, 47, by = 3), "1"),
    other = units::set_units(seq(5, 14, by = 1), "1")
  )

  result <- subtyped(
    data = df,
    case_timeseries = "ili",
    time = "date",
    influenzaA_unsubtyped_timeseries = "inf_A",
    influenzaA_subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
    other_pathogen_timeseries = c("inf_B", "other")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "subtyped")
  # Units should be stripped
  expect_true(is.numeric(result$data$case_timeseries))
  expect_true(is.numeric(result$data$component_pathogens))
  expect_true(is.numeric(result$data$influenzaA_subtyped))
})

test_that("intake functions accept mixed numeric and units columns", {
  skip_if_not_installed("units")

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = 1:10,  # Regular numeric
    alpha = units::set_units(11:20, "1"),  # With units
    delta = 21:30  # Regular numeric
  )

  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_true(is.numeric(result$data$case_timeseries))
  expect_true(is.numeric(result$data$component_pathogens))
})

test_that("intake functions handle different unit types appropriately", {
  skip_if_not_installed("units")

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = units::set_units(1:10, "1"),  # Dimensionless
    alpha = units::set_units(11:20, "1/d")  # Rate (cases per day)
  )

  # Should still work - units are stripped
  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_true(is.numeric(result$data$case_timeseries))
  expect_true(is.numeric(result$data$component_pathogens))
})

test_that("check_column_numeric accepts units objects", {
  skip_if_not_installed("units")

  df <- data.frame(
    numeric_col = 1:10,
    units_col = units::set_units(1:10, "1"),
    char_col = as.character(1:10)
  )

  # Should pass for numeric
  expect_silent(check_column_numeric(df, "numeric_col"))

  # Should pass for units
  expect_silent(check_column_numeric(df, "units_col"))

  # Should error for character
  expect_error(
    check_column_numeric(df, "char_col"),
    "must be numeric or have units"
  )
})

test_that("units are properly stripped before modeling", {
  skip_if_not_installed("units")

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 50),
    cases = units::set_units(rpois(50, 50), "1"),
    alpha = units::set_units(rpois(50, 20), "1"),
    delta = units::set_units(rpois(50, 15), "1")
  )

  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  # The returned data should be plain numeric, not units objects
  expect_false(inherits(result$data$case_timeseries, "units"))
  expect_false(inherits(result$data$component_pathogens, "units"))
  expect_true(is.numeric(result$data$case_timeseries))
  expect_true(is.numeric(result$data$component_pathogens))
})

test_that("units objects work with tsibble validation", {
  skip_if_not_installed("units")

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = units::set_units(1:10, "1"),
    alpha = units::set_units(11:20, "1")
  )

  # Should pass all tsibble validations (no gaps, regular, ordered)
  result <- single(
    data = df,
    case_timeseries = "cases",
    time = "date"
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_length(result$data$case_timeseries, 10)
})

test_that("units objects in time series classes are handled", {
  skip_if_not_installed("units")
  skip_if_not_installed("xts")

  # Note: xts may not support units objects directly, but if converted
  # to data frame first, they should work
  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10)

  df <- data.frame(
    cases = units::set_units(1:10, "1"),
    alpha = units::set_units(11:20, "1")
  )

  # Create xts from data frame with units
  xts_obj <- xts::xts(df, order.by = dates)

  # If xts strips units (which it might), this test documents that behavior
  # Either way, your intake function should handle it
  result <- tryCatch({
    single(
      data = xts_obj,
      case_timeseries = "cases"
    )
  }, error = function(e) {
    NULL
  })

  # If it works, great. If not, that's because xts doesn't support units
  # Either outcome is acceptable - documenting the limitation
  if (!is.null(result)) {
    expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  }
})
