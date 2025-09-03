# Pathogen structure tests using helper functions

# Tests for single() function
test_that("single() creates correct structure with default and custom pathogen names", {
  check_package_data()

  # Test with default pathogen name
  result_default <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date
  )

  expect_s3_class(result_default, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result_default$pathogen_structure, "single")
  expect_equal(result_default$pathogen_names, "default")
  expect_equal(result_default$data$case_timeseries, sarscov2$cases)
  expect_equal(result_default$data$time, sarscov2$date)

  # Test with custom pathogen name
  result_custom <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "SARS-CoV-2"
  )

  expect_equal(result_custom$pathogen_names, "SARS-CoV-2")
  expect_equal(result_custom$data$case_timeseries, sarscov2$cases)
})

# Tests for multiple() function
test_that("multiple() creates correct structure with different parameter combinations", {
  check_package_data()
  expected_names <- get_expected_pathogen_names()

  test_cases <- list(
    list(
      name = "default_parameters",
      smoothing = "shared", noise = "observation_noise_only",
      expected_cov = 0, expected_noise = 0
    ),
    list(
      name = "independent_smoothing",
      smoothing = "independent", noise = "observation_noise_only",
      expected_cov = 1, expected_noise = 0
    ),
    list(
      name = "correlated_smoothing",
      smoothing = "correlated", noise = "pathogen_specific_noise",
      expected_cov = 2, expected_noise = 1
    )
  )

  for (test_case in test_cases) {
    result <- multiple(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date,
      component_pathogen_timeseries = list(
        alpha = sarscov2$alpha,
        delta = sarscov2$delta,
        omicron = sarscov2$omicron,
        other = sarscov2$other
      ),
      smoothing_structure = test_case$smoothing,
      observation_noise = test_case$noise
    )

    # Test basic structure
    expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
    expect_equal(result$pathogen_structure, "multiple")
    expect_equal(result$pathogen_names, expected_names$sarscov2_multiple)

    # Test data structure
    expect_equal(result$data$case_timeseries, sarscov2$cases)
    expect_equal(result$data$time, sarscov2$date)
    expect_equal(nrow(result$data$component_pathogens), 4)
    expect_equal(ncol(result$data$component_pathogens), length(sarscov2$cases))

    # Test model parameters
    expect_equal(result$model_params$cov_structure, test_case$expected_cov,
                 info = paste("Failed for", test_case$name))
    expect_equal(result$model_params$noise_structure, test_case$expected_noise,
                 info = paste("Failed for", test_case$name))
  }
})

test_that("multiple() correctly handles component pathogen data", {
  check_package_data()

  # Test matrix transposition and single pathogen case
  result_multi <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta
    )
  )

  # Check matrix structure and transposition
  expect_equal(result_multi$data$component_pathogens[1, ], sarscov2$alpha)
  expect_equal(result_multi$data$component_pathogens[2, ], sarscov2$delta)

  # Test single pathogen in component list
  result_single <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(alpha = sarscov2$alpha)
  )

  expect_equal(result_single$pathogen_names, "alpha")
  expect_equal(nrow(result_single$data$component_pathogens), 1)
  expect_equal(ncol(result_single$data$component_pathogens), length(sarscov2$cases))
})

# Tests for subtyped() function
test_that("subtyped() creates correct structure with different parameter combinations", {
  check_package_data()
  expected_names <- get_expected_pathogen_names()

  test_cases <- list(
    list(
      name = "default_parameters",
      smoothing = "shared", noise = "observation_noise_only",
      expected_cov = 0, expected_noise = 0
    ),
    list(
      name = "independent_smoothing",
      smoothing = "independent", noise = "pathogen_specific_noise",
      expected_cov = 1, expected_noise = 1
    )
  )

  for (test_case in test_cases) {
    result <- subtyped(
      case_timeseries = influenza$ili,
      time = influenza$week,
      influenzaA_unsubtyped_timeseries = influenza$inf_A,
      influenzaA_subtyped_timeseries = list(
        H3N2 = influenza$inf_H3N2,
        H1N1 = influenza$inf_H1N1
      ),
      other_pathogen_timeseries = list(
        influenzaB = influenza$inf_B
      ),
      smoothing_structure = test_case$smoothing,
      observation_noise = test_case$noise
    )

    # Test basic structure
    expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
    expect_equal(result$pathogen_structure, "subtyped")
    expect_equal(result$pathogen_names, c(expected_names$influenza_subtyped, "influenzaB"))

    # Test data structure
    expect_equal(result$data$case_timeseries, influenza$ili)
    expect_equal(result$data$time, influenza$week)
    expect_equal(nrow(result$data$component_pathogens), 2) # influenzaA + influenzaB
    expect_equal(ncol(result$data$component_pathogens), length(influenza$ili))
    expect_equal(nrow(result$data$influenzaA_subtyped), 2) # H3N2 + H1N1
    expect_equal(ncol(result$data$influenzaA_subtyped), length(influenza$ili))

    # Test model parameters
    expect_equal(result$model_params$cov_structure, test_case$expected_cov,
                 info = paste("Failed for", test_case$name))
    expect_equal(result$model_params$noise_structure, test_case$expected_noise,
                 info = paste("Failed for", test_case$name))
  }
})

test_that("subtyped() correctly handles matrix transposition and pathogen ordering", {
  check_package_data()
  expected_names <- get_expected_pathogen_names()

  result <- subtyped(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(
      H3N2 = influenza$inf_H3N2,
      H1N1 = influenza$inf_H1N1
    ),
    other_pathogen_timeseries = list(
      influenzaB = influenza$inf_B,
      other = influenza$num_spec - influenza$inf_all
    )
  )

  # Test pathogen name ordering
  expected_full_names <- c(expected_names$influenza_subtyped, expected_names$influenza_other)
  expect_equal(result$pathogen_names, expected_full_names)

  # Test matrix transposition
  expect_equal(result$data$component_pathogens[1, ], influenza$inf_A) # influenzaA
  expect_equal(result$data$component_pathogens[2, ], influenza$inf_B) # influenzaB
  expect_equal(result$data$influenzaA_subtyped[1, ], influenza$inf_H3N2) # H3N2
  expect_equal(result$data$influenzaA_subtyped[2, ], influenza$inf_H1N1) # H1N1
})

# Tests for utility functions - these are simple and can stay as-is
test_that("utility functions return correct mappings", {
  # Test get_cov_structure()
  expect_equal(get_cov_structure('shared'), 0)
  expect_equal(get_cov_structure('independent'), 1)
  expect_equal(get_cov_structure('correlated'), 2)

  # Test get_noise_structure()
  expect_equal(get_noise_structure('observation_noise_only'), 0)
  expect_equal(get_noise_structure('pathogen_specific_noise'), 1)
})

# Simplified validation tests using helper function
validate_length_mismatch <- function(constructor_func, base_args, mismatch_args, expected_pattern = "length|same|match") {
  args <- modifyList(base_args, mismatch_args)
  expect_error(do.call(constructor_func, args), expected_pattern, ignore.case = TRUE)
}

test_that("all pathogen structure functions validate input vector lengths", {
  check_package_data()

  # Test single() length validation
  single_base <- list(case_timeseries = sarscov2$cases, time = sarscov2$date)

  validate_length_mismatch(single, single_base,
                           list(time = sarscov2$date[1:5]))
  validate_length_mismatch(single, single_base,
                           list(time = c(sarscov2$date, sarscov2$date[1:3])))

  # Test multiple() length validation
  multiple_base <- list(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(alpha = sarscov2$alpha)
  )

  validate_length_mismatch(multiple, multiple_base,
                           list(time = sarscov2$date[1:5]))
  validate_length_mismatch(multiple, multiple_base,
                           list(component_pathogen_timeseries = list(alpha = sarscov2$alpha[1:5])))
  validate_length_mismatch(multiple, multiple_base,
                           list(component_pathogen_timeseries = list(
                             alpha = sarscov2$alpha,
                             delta = sarscov2$delta[1:5]
                           )))

  # Test subtyped() length validation
  subtyped_base <- list(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
    other_pathogen_timeseries = list(B = influenza$inf_B)
  )

  validate_length_mismatch(subtyped, subtyped_base,
                           list(time = influenza$week[1:10]))
  validate_length_mismatch(subtyped, subtyped_base,
                           list(influenzaA_unsubtyped_timeseries = influenza$inf_A[1:10]))
  validate_length_mismatch(subtyped, subtyped_base,
                           list(influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2[1:10])))
  validate_length_mismatch(subtyped, subtyped_base,
                           list(other_pathogen_timeseries = list(B = influenza$inf_B[1:10])))
})

test_that("pathogen structure functions validate parameter arguments", {
  check_package_data()

  # Test multiple() parameter validation
  base_multiple_args <- list(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(alpha = sarscov2$alpha)
  )

  expect_error(do.call(multiple, modifyList(base_multiple_args,
                                            list(smoothing_structure = 'invalid_option'))))
  expect_error(do.call(multiple, modifyList(base_multiple_args,
                                            list(observation_noise = 'invalid_option'))))

  # Test subtyped() parameter validation
  base_subtyped_args <- list(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
    other_pathogen_timeseries = list(B = influenza$inf_B)
  )

  expect_error(do.call(subtyped, modifyList(base_subtyped_args,
                                            list(smoothing_structure = 'invalid_option'))))
  expect_error(do.call(subtyped, modifyList(base_subtyped_args,
                                            list(observation_noise = 'invalid_option'))))
})

test_that("pathogen structure functions preserve data types and handle edge cases", {
  check_package_data()

  # Test data type preservation
  int_cases <- as.integer(sarscov2$cases)
  result <- single(case_timeseries = int_cases, time = sarscov2$date)
  expect_type(result$data$case_timeseries, "integer")
  expect_s3_class(result$data$time, "Date")

  # Test pathogen name preservation from named lists
  pathogen_names <- c("variant_alpha", "variant_delta", "variant_omicron")
  component_list <- list(
    variant_alpha = sarscov2$alpha,
    variant_delta = sarscov2$delta,
    variant_omicron = sarscov2$omicron
  )

  result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = component_list
  )

  expect_equal(result$pathogen_names, pathogen_names)
})

test_that("all pathogen structure functions return required components", {
  check_package_data()

  # Test single() components
  single_result <- single(case_timeseries = sarscov2$cases, time = sarscov2$date)
  expect_true(all(c("pathogen_structure", "pathogen_names", "data") %in% names(single_result)))
  expect_true(all(c("case_timeseries", "time") %in% names(single_result$data)))

  # Test multiple() components
  multiple_result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(alpha = sarscov2$alpha)
  )
  expect_true(all(c("pathogen_structure", "pathogen_names", "data", "model_params") %in% names(multiple_result)))
  expect_true(all(c("case_timeseries", "time", "component_pathogens") %in% names(multiple_result$data)))
  expect_true(all(c("cov_structure", "noise_structure") %in% names(multiple_result$model_params)))

  # Test subtyped() components
  subtyped_result <- subtyped(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
    other_pathogen_timeseries = list(B = influenza$inf_B)
  )
  expect_true(all(c("pathogen_structure", "pathogen_names", "data", "model_params") %in% names(subtyped_result)))
  expect_true(all(c("case_timeseries", "time", "component_pathogens", "influenzaA_subtyped") %in% names(subtyped_result$data)))
})
