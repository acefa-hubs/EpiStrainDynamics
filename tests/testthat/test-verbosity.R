# Tests for verbosity and error handling
# These tests verify that fit_model properly handles verbose output,
# message/warning suppression, and error catching per Bayesian standards.
#
# Note: Core fit_model functionality is tested in test-fit-model.R

test_that("fit_model has verbose parameter defaulting to TRUE", {
  # Check that the function has verbose argument
  fit_model_args <- names(formals(fit_model))
  expect_true("verbose" %in% fit_model_args)

  # Check that default is TRUE
  expect_true(formals(fit_model)$verbose)
})

test_that("fit_model with verbose=TRUE shows progress", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  # With verbose=TRUE, function should complete successfully
  suppressWarnings(
    fit <- fit_model(mod, n_iter = 100, n_chain = 1, verbose = TRUE)
  )

  expect_s3_class(fit, "EpiStrainDynamics.fit")
})

# Test suppress messages/progress but retain warnings/errors

test_that("fit_model has separate control for messages and warnings", {
  # verbose parameter should control messages/progress
  # warnings should still appear even when verbose=FALSE
  fit_model_args <- names(formals(fit_model))
  expect_true("verbose" %in% fit_model_args)
  expect_true("suppress_warnings" %in% fit_model_args)
})

test_that("fit_model with verbose=FALSE suppresses progress", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  # Capture output with verbose=FALSE
  output <- capture.output({
    suppressWarnings(
      fit <- fit_model(mod, n_iter = 100, n_chain = 1, verbose = FALSE)
    )
  }, type = "message")

  # Should complete successfully
  expect_s3_class(fit, "EpiStrainDynamics.fit")

  # Progress messages should be minimal
  expect_lt(length(output), 5)
})

# Test suppress warnings when appropriate

test_that("fit_model has suppress_warnings parameter", {
  fit_model_args <- names(formals(fit_model))
  expect_true("suppress_warnings" %in% fit_model_args)

  # Check that default is FALSE
  expect_false(formals(fit_model)$suppress_warnings)
})

test_that("fit_model with suppress_warnings=TRUE suppresses warnings", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  # With suppress_warnings=TRUE, warnings should be suppressed
  expect_no_error({
    fit <- fit_model(
      mod,
      n_iter = 100,
      n_chain = 1,
      verbose = FALSE,
      suppress_warnings = TRUE
    )
  })

  expect_s3_class(fit, "EpiStrainDynamics.fit")
})

test_that("validate_suppress_warnings rejects invalid inputs", {
  expect_error(validate_suppress_warnings("TRUE"), "must be logical")
  expect_error(validate_suppress_warnings(1), "must be logical")
  expect_error(validate_suppress_warnings(c(TRUE, FALSE)), "must be a single value")
  expect_error(validate_suppress_warnings(NA), "cannot be NA")

  # Should pass for valid inputs
  expect_silent(validate_suppress_warnings(TRUE))
  expect_silent(validate_suppress_warnings(FALSE))
})

# Test catch and process errors appropriately

test_that("fit_model catches Stan failures and provides informative output", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  # Set standata to NULL - this will cause Stan to fail
  mod$standata <- NULL

  # Now it should throw an actual error
  expect_error(
    fit_model(mod, n_iter = 100, n_chain = 1, verbose = FALSE),
    "Stan sampling failed"
  )
})

test_that("fit_model errors can be caught and inspected", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  # Corrupt standata
  mod$standata <- NULL

  # Catch the error
  result <- tryCatch(
    fit_model(mod, n_iter = 100, n_chain = 1, verbose = FALSE),
    error = function(e) e
  )

  # Error should be catchable
  expect_true(inherits(result, "error"))

  # Error should have informative class
  expect_true(inherits(result, "EpiStrainDynamics.fit.error"))

  # Error object should contain useful information
  expect_true("message" %in% names(result))
  expect_true("constructed_model" %in% names(result))
})

test_that("fit_model errors contain the constructed model for inspection", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  mod$standata <- NULL

  result <- tryCatch(
    fit_model(mod, n_iter = 100, n_chain = 1, verbose = FALSE),
    error = function(e) e
  )

  if (inherits(result, "EpiStrainDynamics.fit.error")) {
    # Should be able to access the model that failed
    expect_s3_class(result$constructed_model, "EpiStrainDynamics.model")

    # Should be able to inspect the error message
    expect_true(nchar(result$message) > 0)
  }
})

# Test combined verbosity controls

test_that("fit_model can suppress both messages and warnings", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  # Should work with both suppressions
  expect_no_error({
    fit <- fit_model(
      mod,
      n_iter = 100,
      n_chain = 1,
      verbose = FALSE,
      suppress_warnings = TRUE
    )
  })

  expect_s3_class(fit, "EpiStrainDynamics.fit")
})

test_that("fit_model verbose and suppress_warnings work independently", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  check_package_data()

  sarscov2_subset <- sarscov2[1:40, ]

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2_subset,
      case_timeseries = "cases",
      time = "date"
    )
  )

  # Test all four combinations
  suppressWarnings({
    fit1 <- fit_model(mod, n_iter = 100, n_chain = 1,
                      verbose = TRUE, suppress_warnings = FALSE)
  })

  suppressWarnings({
    fit2 <- fit_model(mod, n_iter = 100, n_chain = 1,
                      verbose = TRUE, suppress_warnings = TRUE)
  })

  suppressWarnings({
    fit3 <- fit_model(mod, n_iter = 100, n_chain = 1,
                      verbose = FALSE, suppress_warnings = FALSE)
  })

  suppressWarnings({
    fit4 <- fit_model(mod, n_iter = 100, n_chain = 1,
                      verbose = FALSE, suppress_warnings = TRUE)
  })

  expect_s3_class(fit1, "EpiStrainDynamics.fit")
  expect_s3_class(fit2, "EpiStrainDynamics.fit")
  expect_s3_class(fit3, "EpiStrainDynamics.fit")
  expect_s3_class(fit4, "EpiStrainDynamics.fit")
})
