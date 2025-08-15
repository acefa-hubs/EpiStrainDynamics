# Load package datasets for testing
data(sarscov2, envir = environment())
data(influenza, envir = environment())

# Verify datasets loaded successfully
if (!exists("sarscov2") || !exists("influenza")) {
  stop("Required datasets (sarscov2, influenza) could not be loaded")
}

# Tests for single() function
test_that("single() creates correct structure with default pathogen name", {
  result <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "single")
  expect_equal(result$pathogen_names, "default")
  expect_equal(result$data$case_timeseries, sarscov2$cases)
  expect_equal(result$data$time, sarscov2$date)
})

test_that("single() creates correct structure with custom pathogen name", {
  result <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    pathogen_name = "SARS-CoV-2"
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "single")
  expect_equal(result$pathogen_names, "SARS-CoV-2")
  expect_equal(result$data$case_timeseries, sarscov2$cases)
  expect_equal(result$data$time, sarscov2$date)
})

# Tests for multiple() function
test_that("multiple() creates correct structure with default parameters", {
  result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta,
      omicron = sarscov2$omicron
    )
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "multiple")
  expect_equal(result$pathogen_names, c("alpha", "delta", "omicron"))
  expect_equal(result$data$case_timeseries, sarscov2$cases)
  expect_equal(result$data$time, sarscov2$date)
  expect_equal(nrow(result$data$component_pathogens), 3) # 3 pathogens
  expect_equal(ncol(result$data$component_pathogens), length(sarscov2$cases)) # same length as cases
  expect_equal(result$model_params$cov_structure, 0) # default 'shared'
  expect_equal(result$model_params$noise_structure, 0) # default 'observation_noise_only'
})

test_that("multiple() works with independent smoothing structure", {
  result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta
    ),
    smoothing_structure = 'independent',
    observation_noise = 'observation_noise_only'
  )

  expect_equal(result$model_params$cov_structure, 1) # independent
  expect_equal(result$model_params$noise_structure, 0) # observation_noise_only
})

test_that("multiple() works with correlated smoothing structure", {
  result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta
    ),
    smoothing_structure = 'correlated',
    observation_noise = 'pathogen_specific_noise'
  )

  expect_equal(result$model_params$cov_structure, 2) # correlated
  expect_equal(result$model_params$noise_structure, 1) # pathogen_specific_noise
})

test_that("multiple() correctly transposes component pathogen matrix", {
  result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha,
      delta = sarscov2$delta
    )
  )

  # Check that first row corresponds to first pathogen (alpha)
  expect_equal(result$data$component_pathogens[1, ], sarscov2$alpha)
  # Check that second row corresponds to second pathogen (delta)
  expect_equal(result$data$component_pathogens[2, ], sarscov2$delta)
})

test_that("multiple() handles single pathogen in component list", {
  result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(
      alpha = sarscov2$alpha
    )
  )

  expect_equal(result$pathogen_names, "alpha")
  expect_equal(nrow(result$data$component_pathogens), 1)
  expect_equal(ncol(result$data$component_pathogens), length(sarscov2$cases))
})

# Tests for subtyped() function
test_that("subtyped() creates correct structure with default parameters", {
  result <- subtyped(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(
      influenzaA.H3N2 = influenza$inf_H3N2,
      influenzaA.H1N1 = influenza$inf_H1N1
    ),
    other_pathogen_timeseries = list(
      influenzaB = influenza$inf_B
    )
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "subtyped")
  expect_equal(result$pathogen_names, c("influenzaA.H3N2", "influenzaA.H1N1", "influenzaB"))
  expect_equal(result$data$case_timeseries, influenza$ili)
  expect_equal(result$data$time, influenza$week)
  expect_equal(nrow(result$data$component_pathogens), 2) # influenzaA + influenzaB
  expect_equal(ncol(result$data$component_pathogens), length(influenza$ili))
  expect_equal(nrow(result$data$influenzaA_subtyped), 2) # H3N2 + H1N1
  expect_equal(ncol(result$data$influenzaA_subtyped), length(influenza$ili))
})

test_that("subtyped() works with independent smoothing structure", {
  result <- subtyped(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(
      influenzaA.H3N2 = influenza$inf_H3N2,
      influenzaA.H1N1 = influenza$inf_H1N1
    ),
    other_pathogen_timeseries = list(
      influenzaB = influenza$inf_B
    ),
    smoothing_structure = 'independent',
    observation_noise = 'pathogen_specific_noise'
  )

  expect_equal(result$model_params$cov_structure, 1) # independent
  expect_equal(result$model_params$noise_structure, 1) # pathogen_specific_noise
})

test_that("subtyped() correctly orders pathogen names", {
  result <- subtyped(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(
      influenzaA.H3N2 = influenza$inf_H3N2,
      influenzaA.H1N1 = influenza$inf_H1N1
    ),
    other_pathogen_timeseries = list(
      influenzaB = influenza$inf_B,
      other = influenza$num_spec - influenza$inf_all
    )
  )

  expected_names <- c("influenzaA.H3N2", "influenzaA.H1N1", "influenzaB", "other")
  expect_equal(result$pathogen_names, expected_names)
})

test_that("subtyped() correctly transposes matrices", {
  result <- subtyped(
    case_timeseries = influenza$ili,
    time = influenza$week,
    influenzaA_unsubtyped_timeseries = influenza$inf_A,
    influenzaA_subtyped_timeseries = list(
      influenzaA.H3N2 = influenza$inf_H3N2,
      influenzaA.H1N1 = influenza$inf_H1N1
    ),
    other_pathogen_timeseries = list(
      influenzaB = influenza$inf_B
    )
  )

  # Check component_pathogens matrix (influenzaA + other pathogens)
  expect_equal(result$data$component_pathogens[1, ], influenza$inf_A) # influenzaA
  expect_equal(result$data$component_pathogens[2, ], influenza$inf_B) # influenzaB

  # Check influenzaA_subtyped matrix
  expect_equal(result$data$influenzaA_subtyped[1, ], influenza$inf_H3N2) # H3N2
  expect_equal(result$data$influenzaA_subtyped[2, ], influenza$inf_H1N1) # H1N1
})

# Tests for get_cov_structure() function
test_that("get_cov_structure() returns correct values", {
  expect_equal(get_cov_structure('shared'), 0)
  expect_equal(get_cov_structure('independent'), 1)
  expect_equal(get_cov_structure('correlated'), 2)
})

# Tests for get_noise_structure() function
test_that("get_noise_structure() returns correct values", {
  expect_equal(get_noise_structure('observation_noise_only'), 0)
  expect_equal(get_noise_structure('pathogen_specific_noise'), 1)
})

# Length validation tests
test_that("single() validates that case_timeseries and time have same length", {
  # Test with shorter time vector - should throw error
  shorter_time <- sarscov2$date[1:5]

  expect_error({
    single(
      case_timeseries = sarscov2$cases,
      time = shorter_time
    )
  }, regexp = "length|same|match", ignore.case = TRUE)

  # Test with longer time vector - should throw error
  longer_time <- c(sarscov2$date, sarscov2$date[1:3])

  expect_error({
    single(
      case_timeseries = sarscov2$cases,
      time = longer_time
    )
  }, regexp = "length|same|match", ignore.case = TRUE)
})

test_that("multiple() validates that all input vectors have same length", {
  # Test case_timeseries vs time mismatch
  shorter_time <- sarscov2$date[1:5]

  expect_error({
    multiple(
      case_timeseries = sarscov2$cases,
      time = shorter_time,
      component_pathogen_timeseries = list(alpha = sarscov2$alpha)
    )
  }, regexp = "length|same|match", ignore.case = TRUE)

  # Test case_timeseries vs component_pathogen_timeseries mismatch
  shorter_alpha <- sarscov2$alpha[1:5]

  expect_error({
    multiple(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date,
      component_pathogen_timeseries = list(alpha = shorter_alpha)
    )
  }, regexp = "length|same|match", ignore.case = TRUE)

  # Test mismatched component pathogen lengths
  shorter_delta <- sarscov2$delta[1:5]

  expect_error({
    multiple(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date,
      component_pathogen_timeseries = list(
        alpha = sarscov2$alpha,
        delta = shorter_delta  # Different length
      )
    )
  }, regexp = "length|same|match", ignore.case = TRUE)
})

test_that("subtyped() validates that all input vectors have same length", {
  # Test case_timeseries vs time mismatch
  shorter_time <- influenza$week[1:10]

  expect_error({
    subtyped(
      case_timeseries = influenza$ili,
      time = shorter_time,
      influenzaA_unsubtyped_timeseries = influenza$inf_A,
      influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
      other_pathogen_timeseries = list(B = influenza$inf_B)
    )
  }, regexp = "length|same|match", ignore.case = TRUE)

  # Test case_timeseries vs influenzaA_unsubtyped_timeseries mismatch
  shorter_inf_A <- influenza$inf_A[1:10]

  expect_error({
    subtyped(
      case_timeseries = influenza$ili,
      time = influenza$week,
      influenzaA_unsubtyped_timeseries = shorter_inf_A,
      influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
      other_pathogen_timeseries = list(B = influenza$inf_B)
    )
  }, regexp = "length|same|match", ignore.case = TRUE)

  # Test mismatched influenzaA_subtyped_timeseries lengths
  shorter_H3N2 <- influenza$inf_H3N2[1:10]

  expect_error({
    subtyped(
      case_timeseries = influenza$ili,
      time = influenza$week,
      influenzaA_unsubtyped_timeseries = influenza$inf_A,
      influenzaA_subtyped_timeseries = list(
        H3N2 = shorter_H3N2,  # Different length
        H1N1 = influenza$inf_H1N1
      ),
      other_pathogen_timeseries = list(B = influenza$inf_B)
    )
  }, regexp = "length|same|match", ignore.case = TRUE)

  # Test mismatched other_pathogen_timeseries lengths
  shorter_inf_B <- influenza$inf_B[1:10]

  expect_error({
    subtyped(
      case_timeseries = influenza$ili,
      time = influenza$week,
      influenzaA_unsubtyped_timeseries = influenza$inf_A,
      influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
      other_pathogen_timeseries = list(B = shorter_inf_B)  # Different length
    )
  }, regexp = "length|same|match", ignore.case = TRUE)

  # Test mismatched lengths within other_pathogen_timeseries
  other_short <- rep(10, 10)  # Different length

  expect_error({
    subtyped(
      case_timeseries = influenza$ili,
      time = influenza$week,
      influenzaA_unsubtyped_timeseries = influenza$inf_A,
      influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
      other_pathogen_timeseries = list(
        B = influenza$inf_B,
        other = other_short  # Different length
      )
    )
  }, regexp = "length|same|match", ignore.case = TRUE)
})

test_that("pathogen structure functions preserve data types", {
  # Test with integer case data
  int_cases <- as.integer(sarscov2$cases)
  result <- single(case_timeseries = int_cases, time = sarscov2$date)
  expect_type(result$data$case_timeseries, "integer")

  # Test with Date time data
  result <- single(case_timeseries = sarscov2$cases, time = sarscov2$date)
  expect_s3_class(result$data$time, "Date")
})

test_that("multiple() preserves pathogen names from named list", {
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

# Test to ensure all required components are present
test_that("all pathogen structure functions return required components", {
  # Test single()
  single_result <- single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date
  )

  expect_true(all(c("pathogen_structure", "pathogen_names", "data") %in% names(single_result)))
  expect_true(all(c("case_timeseries", "time") %in% names(single_result$data)))

  # Test multiple()
  multiple_result <- multiple(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date,
    component_pathogen_timeseries = list(alpha = sarscov2$alpha)
  )

  expect_true(all(c("pathogen_structure", "pathogen_names", "data", "model_params") %in% names(multiple_result)))
  expect_true(all(c("case_timeseries", "time", "component_pathogens") %in% names(multiple_result$data)))
  expect_true(all(c("cov_structure", "noise_structure") %in% names(multiple_result$model_params)))

  # Test subtyped()
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

# Test argument matching for smoothing_structure and observation_noise
test_that("multiple() and subtyped() validate smoothing_structure arguments", {
  # Test that invalid smoothing_structure throws error
  expect_error({
    multiple(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date,
      component_pathogen_timeseries = list(alpha = sarscov2$alpha),
      smoothing_structure = 'invalid_option'
    )
  })

  # Test that invalid observation_noise throws error
  expect_error({
    multiple(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date,
      component_pathogen_timeseries = list(alpha = sarscov2$alpha),
      observation_noise = 'invalid_option'
    )
  })
})

test_that("subtyped() validates smoothing_structure and observation_noise arguments", {
  # Test that invalid smoothing_structure throws error
  expect_error({
    subtyped(
      case_timeseries = influenza$ili,
      time = influenza$week,
      influenzaA_unsubtyped_timeseries = influenza$inf_A,
      influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
      other_pathogen_timeseries = list(B = influenza$inf_B),
      smoothing_structure = 'invalid_option'
    )
  })

  # Test that invalid observation_noise throws error
  expect_error({
    subtyped(
      case_timeseries = influenza$ili,
      time = influenza$week,
      influenzaA_unsubtyped_timeseries = influenza$inf_A,
      influenzaA_subtyped_timeseries = list(H3N2 = influenza$inf_H3N2),
      other_pathogen_timeseries = list(B = influenza$inf_B),
      observation_noise = 'invalid_option'
    )
  })
})
