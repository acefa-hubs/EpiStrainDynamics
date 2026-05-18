#' @srrstats {G5.2, G5.2b} error and warning behaviour explicitly tested,
#'   including conditions that trigger these messages and comparing with
#'   expected results
#' @srrstats {G3.0} comparisons made between appropriate values

# ==============================================================================
# TESTS: smoothing_structure() constructor
# ==============================================================================

test_that("smoothing_structure() returns correct class", {
  result <- smoothing_structure("shared")
  expect_s3_class(result, "EpiStrainDynamics.smoothing")
})

test_that("smoothing_structure() accepts all valid types", {
  expect_silent(smoothing_structure("shared"))
  expect_silent(smoothing_structure("independent"))
  expect_silent(smoothing_structure("correlated"))
})

test_that("smoothing_structure() is case insensitive", {
  expect_silent(smoothing_structure("SHARED"))
  expect_silent(smoothing_structure("Independent"))
  expect_silent(smoothing_structure("CORRELATED"))
})

test_that("smoothing_structure() rejects invalid type", {
  expect_error(smoothing_structure("invalid"), class = "rlang_error")
})

test_that("smoothing_structure() alerts when tau priors provided for correlated type", {
  expect_message(
    smoothing_structure("correlated", tau_mean = 0.5, tau_sd = 1.0),
    "ignored for"
  )
})

test_that("smoothing_structure() accepts priors for shared type", {
  result <- smoothing_structure("shared", tau_mean = 0.5, tau_sd = 1.0)
  expect_s3_class(result, "EpiStrainDynamics.smoothing")
  expect_equal(result$smoothing_type, "shared")
  expect_equal(result$priors_provided, 2)
})

test_that("smoothing_structure() accepts priors for independent type", {
  result <- smoothing_structure("independent",
                                tau_mean = c(0.1, 0.2),
                                tau_sd = c(1.0, 1.0))
  expect_s3_class(result, "EpiStrainDynamics.smoothing")
  expect_equal(result$smoothing_type, "independent")
  expect_equal(result$priors_provided, 2)
})

test_that("smoothing_structure() sets priors_provided = 1 when no priors given", {
  result <- smoothing_structure("shared")
  expect_equal(result$priors_provided, 1)
})

# ==============================================================================
# TESTS: validate_smoothing_structure()
# ==============================================================================

test_that("validate_smoothing_structure() errors on wrong class", {
  bad_obj <- list(smoothing_type = "shared", tau_priors = list(), priors_provided = 1)
  expect_error(
    EpiStrainDynamics:::validate_smoothing_structure(bad_obj, "single"),
    "smoothing_params must be created using the smoothing_structure\\(\\) function"
  )
})

test_that("validate_smoothing_structure() returns smoothing object", {
  obj <- smoothing_structure("shared")
  result <- EpiStrainDynamics:::validate_smoothing_structure(obj, "single")
  expect_s3_class(result, "EpiStrainDynamics.smoothing")
})

# ==============================================================================
# TESTS: validate_smoothing_for_single()
# ==============================================================================

test_that("validate_smoothing_for_single() passes through non-single pathogen", {
  obj <- smoothing_structure("independent")
  result <- EpiStrainDynamics:::validate_smoothing_for_single(obj, "multiple")
  expect_equal(result$smoothing_type, "independent")
})

test_that("validate_smoothing_for_single() coerces non-shared to shared for single", {
  obj <- smoothing_structure("independent")
  expect_message(
    result <- EpiStrainDynamics:::validate_smoothing_for_single(obj, "single"),
    'smoothing_type can only be "shared"'
  )
  expect_equal(result$smoothing_type, "shared")
})

test_that("validate_smoothing_for_single() sets default priors for single", {
  obj <- smoothing_structure("shared")
  result <- EpiStrainDynamics:::validate_smoothing_for_single(obj, "single")
  expect_equal(result$tau_priors$mean, 0.0)
  expect_equal(result$tau_priors$sd, 1.0)
})

test_that("validate_smoothing_for_single() preserves existing priors for single", {
  obj <- smoothing_structure("shared", tau_mean = 0.5, tau_sd = 2.0)
  result <- EpiStrainDynamics:::validate_smoothing_for_single(obj, "single")
  expect_equal(result$tau_priors$mean, 0.5)
  expect_equal(result$tau_priors$sd, 2.0)
})

# ==============================================================================
# TESTS: validate_smoothing_independent()
# ==============================================================================

test_that("validate_smoothing_independent() passes through non-independent type", {
  obj <- smoothing_structure("shared")
  result <- EpiStrainDynamics:::validate_smoothing_independent(obj, NULL)
  expect_equal(result$smoothing_type, "shared")
})

test_that("validate_smoothing_independent() errors when pathogen_names missing", {
  obj <- smoothing_structure("independent")
  expect_error(
    EpiStrainDynamics:::validate_smoothing_independent(obj, NULL),
    "pathogen_names is required"
  )
})

test_that("validate_smoothing_independent() sets default priors when none provided", {
  obj <- smoothing_structure("independent")
  result <- EpiStrainDynamics:::validate_smoothing_independent(
    obj, c("alpha", "delta", "omicron", "other"))
  expect_equal(result$tau_priors$mean, rep(0.0, 4))
  expect_equal(result$tau_priors$sd, rep(1.0, 4))
})

test_that("validate_smoothing_independent() expands scalar priors to vector", {
  obj <- smoothing_structure("independent", tau_mean = 0.5, tau_sd = 1.0)
  result <- EpiStrainDynamics:::validate_smoothing_independent(
    obj, c("alpha", "delta", "omicron", "other"))
  expect_length(result$tau_priors$mean, 4)
  expect_true(all(result$tau_priors$mean == 0.5))
})

test_that("validate_smoothing_independent() errors on dimension mismatch", {
  obj <- smoothing_structure("independent",
                             tau_mean = c(0.1, 0.2),
                             tau_sd = c(1.0, 1.0))
  expect_error(
    EpiStrainDynamics:::validate_smoothing_independent(
      obj, c("alpha", "delta", "omicron", "other")),
    "does not match number of pathogens"
  )
})

# ==============================================================================
# TESTS: validate_smoothing_dim()
# ==============================================================================

test_that("validate_smoothing_dim() passes through matching dimensions", {
  result <- EpiStrainDynamics:::validate_smoothing_dim(c(0.1, 0.2, 0.3, 0.4), 4, "tau_mean")
  expect_equal(result, c(0.1, 0.2, 0.3, 0.4))
})

test_that("validate_smoothing_dim() expands scalar to vector", {
  result <- EpiStrainDynamics:::validate_smoothing_dim(0.5, 4, "tau_mean")
  expect_equal(result, rep(0.5, 4))
})

test_that("validate_smoothing_dim() errors on wrong vector length", {
  expect_error(
    EpiStrainDynamics:::validate_smoothing_dim(c(0.1, 0.2), 4, "tau_mean"),
    "does not match number of pathogens"
  )
})

# ==============================================================================
# TESTS: validate_smoothing_shared()
# ==============================================================================

test_that("validate_smoothing_shared() passes through non-shared type", {
  obj <- smoothing_structure("independent")
  result <- EpiStrainDynamics:::validate_smoothing_shared(obj)
  expect_equal(result$smoothing_type, "independent")
})

test_that("validate_smoothing_shared() sets default priors when none provided", {
  obj <- smoothing_structure("shared")
  result <- EpiStrainDynamics:::validate_smoothing_shared(obj)
  expect_equal(as.numeric(result$tau_priors$mean), 0.0)
  expect_equal(as.numeric(result$tau_priors$sd), 1.0)
})

test_that("validate_smoothing_shared() preserves existing priors", {
  obj <- smoothing_structure("shared", tau_mean = 0.5, tau_sd = 2.0)
  result <- EpiStrainDynamics:::validate_smoothing_shared(obj)
  expect_equal(result$tau_priors$mean, 0.5)
  expect_equal(result$tau_priors$sd, 2.0)
})

# ==============================================================================
# TESTS: validate_smoothing_correlated()
# ==============================================================================

test_that("validate_smoothing_correlated() passes through non-correlated type", {
  obj <- smoothing_structure("shared")
  result <- EpiStrainDynamics:::validate_smoothing_correlated(obj)
  expect_equal(result$smoothing_type, "shared")
})

test_that("validate_smoothing_correlated() sets empty priors when none provided", {
  obj <- smoothing_structure("correlated")
  result <- EpiStrainDynamics:::validate_smoothing_correlated(obj)
  expect_equal(result$tau_priors$mean, numeric(0))
  expect_equal(result$tau_priors$sd, numeric(0))
})

test_that("validate_smoothing_correlated() preserves existing priors_provided = 2", {
  obj <- smoothing_structure("correlated")
  obj$priors_provided <- 2
  obj$tau_priors$mean <- c(0.1, 0.2)
  obj$tau_priors$sd <- c(1.0, 1.0)
  result <- EpiStrainDynamics:::validate_smoothing_correlated(obj)
  expect_equal(result$tau_priors$mean, c(0.1, 0.2))
})

