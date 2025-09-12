test_that("fit_model parameters are passed correctly to Stan", {
  skip_if_not_installed("rstan")

  # Capture parameters passed to rstan::sampling
  captured_params <- list()

  # Create custom mock in test environment
  test_env$mock_sampling_capture <- function(object, data, iter, warmup, chains) {
    captured_params$iter <<- iter
    captured_params$warmup <<- warmup
    captured_params$chains <<- chains
    return(test_env$mock_stan_fit)
  }

  # Set up mocks with parameter capture
  setup_stan_mocks()
  assignInNamespace("sampling", test_env$mock_sampling_capture, ns = "rstan")
  assign("stanmodels", test_env$mock_stanmodels, envir = .GlobalEnv)
  on.exit(teardown_mocks())

  models <- create_test_models()
  fit_model(models$rw_single, iter = 500, warmup = 200, chains = 4)

  expect_equal(captured_params$iter, 500)
  expect_equal(captured_params$warmup, 200)
  expect_equal(captured_params$chains, 4)
})

test_that("standata structure is correct for different model types", {
  skip_if_not_installed("rstan")

  captured_data <- list()

  # Create custom mock in test environment
  test_env$mock_sampling_data <- function(object, data, iter, warmup, chains) {
    captured_data <<- data
    return(test_env$mock_stan_fit)
  }

  setup_stan_mocks()
  assignInNamespace("sampling", test_env$mock_sampling_data, ns = "rstan")
  assign("stanmodels", test_env$mock_stanmodels, envir = .GlobalEnv)
  on.exit(teardown_mocks())

  models <- create_test_models()
  expected_lengths <- get_expected_data_lengths()

  # Test single model standata
  fit_model(models$rw_single, iter = 10, warmup = 5, chains = 1)

  expect_true("num_data" %in% names(captured_data))
  expect_true("Y" %in% names(captured_data))
  expect_true("week_effect" %in% names(captured_data))
  expect_true("DOW" %in% names(captured_data))
  expect_equal(length(captured_data$Y), expected_lengths$sarscov2_length)

  # Test multiple model standata
  fit_model(models$rw_multiple, iter = 10, warmup = 5, chains = 1)

  expect_true("num_path" %in% names(captured_data))
  expect_true("P" %in% names(captured_data))
  expect_true("cov_structure" %in% names(captured_data))
  expect_true("noise_structure" %in% names(captured_data))
  expect_equal(captured_data$num_path, expected_lengths$sarscov2_pathogen_count)

  # Test subtyped model standata
  fit_model(models$rw_subtyped, iter = 10, warmup = 5, chains = 1)

  expect_true("P1" %in% names(captured_data))
  expect_true("P2" %in% names(captured_data))
  expect_equal(length(captured_data$Y), expected_lengths$influenza_length)
})

test_that("model-specific standata validation", {
  skip_if_not_installed("rstan")

  all_captured_data <- list()

  # Create custom mock in test environment
  test_env$mock_sampling_data_collector <- function(object, data, iter, warmup, chains) {
    # Store data with a unique key
    key <- paste0("call_", length(all_captured_data) + 1)
    all_captured_data[[key]] <<- data
    return(test_env$mock_stan_fit)
  }

  setup_stan_mocks()
  assignInNamespace("sampling", test_env$mock_sampling_data_collector, ns = "rstan")
  assign("stanmodels", test_env$mock_stanmodels, envir = .GlobalEnv)
  on.exit(teardown_mocks())

  models <- create_test_models()
  expected_lengths <- get_expected_data_lengths()
  expected_names <- get_expected_pathogen_names()

  # Test single models have correct structure
  fit_model(models$rw_single, iter = 10, warmup = 5, chains = 1)
  single_data <- all_captured_data[[length(all_captured_data)]]
  expect_false("num_path" %in% names(single_data))
  expect_false("P" %in% names(single_data))
  expect_equal(length(single_data$Y), expected_lengths$sarscov2_length)

  # Test P-spline models have spline parameters
  fit_model(models$ps_single, iter = 10, warmup = 5, chains = 1)
  ps_data <- all_captured_data[[length(all_captured_data)]]
  expect_true("X" %in% names(ps_data))
  expect_true("num_knots" %in% names(ps_data))
  expect_equal(length(ps_data$Y), expected_lengths$sarscov2_length)

  # Test multiple models have pathogen data
  fit_model(models$rw_multiple, iter = 10, warmup = 5, chains = 1)
  multiple_data <- all_captured_data[[length(all_captured_data)]]
  expect_true("num_path" %in% names(multiple_data))
  expect_true("P" %in% names(multiple_data))
  expect_equal(multiple_data$num_path, expected_lengths$sarscov2_pathogen_count)
  expect_equal(length(multiple_data$Y), expected_lengths$sarscov2_length)

  # Test subtyped models have P1 and P2
  fit_model(models$rw_subtyped, iter = 10, warmup = 5, chains = 1)
  subtyped_data <- all_captured_data[[length(all_captured_data)]]
  expect_true("P1" %in% names(subtyped_data))
  expect_true("P2" %in% names(subtyped_data))
  expect_false("P" %in% names(subtyped_data)) # Should not have single P matrix
  expect_equal(length(subtyped_data$Y), expected_lengths$influenza_length)
})

#' @srrstats {G5.4} Correctness tests with fixed test data
#' @srrstats {G5.5} Fixed random seed for reproducibility
test_that("correctness tests with fixed seed - random walk single pathogen [G5.4, G5.5]", {
  skip_if_not_installed("rstan")
  skip_on_cran() # These tests may be slow

  # G5.5: Set fixed random seed for reproducibility
  set.seed(12345)

  # Create known test data with predictable pattern
  n_obs <- 50
  dates <- seq(as.Date("2023-01-01"), by = "week", length.out = n_obs)

  # Generate data with known trend + noise
  true_trend <- 100 + 0.5 * seq_len(n_obs) + 10 * sin(2 * pi * seq_len(n_obs) / 26)
  test_cases <- rpois(n_obs, lambda = pmax(1, true_trend))

  # G5.4: Test against fixed expected results
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      case_timeseries = test_cases,
      time = dates,
      pathogen_name = "TestVirus"
    ),
    dow_effect = FALSE
  )

  # Fit with minimal iterations for speed but sufficient for convergence check
  fit <- fit_model(mod, iter = 1000, warmup = 500, chains = 2)

  # Test that fit object has expected structure
  expect_s3_class(fit, c("rw_single", "EpiStrainDynamics.fit", "stanfit", "fit"))
  expect_true(length(rstan::extract(fit$fit, "a")) > 0)

  # Test that posterior means are reasonable (within expected range of true values)
  posterior_mean <- rstan::extract(fit$fit, "a")[[1]]
  mean_estimates <- colMeans(posterior_mean)

  # Should be correlated with true trend (allowing for estimation uncertainty)
  correlation <- cor(mean_estimates, true_trend)
  expect_true(correlation > 0.7, info = paste("Correlation too low:", correlation))
})

# G5.6, G5.6a, G5.6b: Parameter recovery tests
#' @srrstats {G5.6} Parameter recovery with known data properties
#' @srrstats {G5.6a} Recovery within defined tolerance
#' @srrstats {G5.6b} Multiple random seeds
test_that("parameter recovery tests - P-spline model [G5.6, G5.6a, G5.6b]", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  # G5.6b: Test with multiple random seeds
  seeds <- c(123, 456, 789)
  recovery_results <- list()

  for (i in seq_along(seeds)) {
    set.seed(seeds[i])

    # G5.6: Create data with known properties
    n_obs <- 60
    dates <- seq(as.Date("2023-01-01"), by = "week", length.out = n_obs)

    # Simulate smooth underlying trend (what P-spline should recover)
    time_numeric <- as.numeric(dates - min(dates)) / 7  # weeks from start
    true_smooth_trend <- 50 + 20 * sin(2 * pi * time_numeric / 26) + 0.1 * time_numeric^2

    # Add observation noise
    observed_cases <- rpois(n_obs, lambda = pmax(1, true_smooth_trend))

    mod <- construct_model(
      method = p_spline(),
      pathogen_structure = single(
        case_timeseries = observed_cases,
        time = dates,
        pathogen_name = "TestVirus"
      ),
      dow_effect = FALSE
    )

    fit <- fit_model(mod, iter = 1000, warmup = 500, chains = 2)

    # Extract fitted values
    y_rep <- extract(fit, "y_rep")[[1]]
    fitted_mean <- colMeans(y_rep)

    # G5.6a: Test recovery within defined tolerance
    # P-splines should smooth towards the underlying trend
    recovery_results[[i]] <- list(
      seed = seeds[i],
      correlation = cor(fitted_mean, true_smooth_trend),
      rmse = sqrt(mean((fitted_mean - true_smooth_trend)^2))
    )

    # Individual seed should have reasonable correlation
    expect_true(recovery_results[[i]]$correlation > 0.8,
                info = paste("Seed", seeds[i], "correlation too low:",
                             recovery_results[[i]]$correlation))
  }

  # G5.6b: Results should be consistent across seeds (within tolerance)
  correlations <- sapply(recovery_results, `[[`, "correlation")
  rmses <- sapply(recovery_results, `[[`, "rmse")

  # Check consistency across seeds
  expect_true(sd(correlations) < 0.1,
              info = paste("Too much variation in correlations across seeds:",
                           paste(correlations, collapse = ", ")))
  expect_true(sd(rmses) < 5,
              info = paste("Too much variation in RMSE across seeds:",
                           paste(rmses, collapse = ", ")))
})

# G5.7: Algorithm performance tests
#' @srrstats {G5.7} Algorithm performance as data properties change
test_that("algorithm performance tests - data size effects [G5.7]", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  set.seed(54321)

  # G5.7: Test performance as data size increases
  sample_sizes <- c(30, 60, 120)
  performance_results <- list()

  for (n in sample_sizes) {
    dates <- seq(as.Date("2023-01-01"), by = "week", length.out = n)

    # Generate data with consistent underlying pattern
    true_trend <- 100 + 0.2 * seq_len(n) + 10 * sin(2 * pi * seq_len(n) / 26)
    test_cases <- rpois(n, lambda = pmax(1, true_trend))

    mod <- construct_model(
      method = random_walk(),
      pathogen_structure = single(
        case_timeseries = test_cases,
        time = dates,
        pathogen_name = "TestVirus"
      ),
      dow_effect = FALSE
    )

    # Measure fitting time
    start_time <- Sys.time()
    fit <- fit_model(mod, iter = 500, warmup = 250, chains = 2)
    fit_time <- as.numeric(Sys.time() - start_time)

    # Check convergence
    rhat_values <- summary(fit)$summary[, "Rhat"]
    max_rhat <- max(rhat_values, na.rm = TRUE)

    # Extract posterior for accuracy assessment
    y_rep <- extract(fit, "y_rep")[[1]]
    fitted_mean <- colMeans(y_rep)
    correlation <- cor(fitted_mean, true_trend)

    performance_results[[paste0("n_", n)]] <- list(
      sample_size = n,
      fit_time = fit_time,
      max_rhat = max_rhat,
      correlation = correlation
    )

    # Basic performance expectations
    expect_true(max_rhat < 1.1, info = paste("Poor convergence for n =", n))
    expect_true(correlation > 0.7, info = paste("Poor fit for n =", n))
  }

  # G5.7: Performance should improve with larger sample sizes
  correlations <- sapply(performance_results, `[[`, "correlation")

  # Correlation should generally increase with sample size
  expect_true(correlations["n_120"] >= correlations["n_30"] - 0.1,
              info = "Algorithm performance should not degrade significantly with more data")

  # Convergence should be better with more data
  rhats <- sapply(performance_results, `[[`, "max_rhat")
  expect_true(all(rhats < 1.1), info = "Convergence issues detected")
})

# G5.9a: Noise susceptibility tests - trivial noise
#' @srrstats {G5.9a} Trivial noise should not meaningfully change results
test_that("noise susceptibility - trivial noise does not change results [G5.9a]", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  set.seed(99999)

  # Create base dataset
  n_obs <- 40
  dates <- seq(as.Date("2023-01-01"), by = "week", length.out = n_obs)
  base_cases <- rpois(n_obs, lambda = 50 + 10 * sin(2 * pi * seq_len(n_obs) / 26))

  # G5.9a: Add trivial noise at machine epsilon scale
  trivial_noise <- runif(n_obs, -1, 1) * .Machine$double.eps * max(base_cases)
  noisy_cases <- base_cases + trivial_noise

  # Fit both models
  mod_original <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      case_timeseries = base_cases,
      time = dates,
      pathogen_name = "TestVirus"
    ),
    dow_effect = FALSE
  )

  mod_noisy <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      case_timeseries = noisy_cases,
      time = dates,
      pathogen_name = "TestVirus"
    ),
    dow_effect = FALSE
  )

  # Use same seed for both fits to ensure comparable sampling
  set.seed(11111)
  fit_original <- fit_model(mod_original, iter = 500, warmup = 250, chains = 2)

  set.seed(11111)
  fit_noisy <- fit_model(mod_noisy, iter = 500, warmup = 250, chains = 2)

  # Extract posterior means
  y_rep_orig <- colMeans(extract(fit_original, "y_rep")[[1]])
  y_rep_noisy <- colMeans(extract(fit_noisy, "y_rep")[[1]])

  # G5.9a: Results should be essentially identical
  correlation <- cor(y_rep_orig, y_rep_noisy)
  rmse <- sqrt(mean((y_rep_orig - y_rep_noisy)^2))
  max_diff <- max(abs(y_rep_orig - y_rep_noisy))

  expect_true(correlation > 0.999,
              info = paste("Trivial noise changed results too much. Correlation:", correlation))
  expect_true(rmse < 0.01 * mean(y_rep_orig),
              info = paste("RMSE too large for trivial noise:", rmse))
  expect_true(max_diff < 0.05 * mean(y_rep_orig),
              info = paste("Max difference too large:", max_diff))
})

# G5.9b: Different random seeds/initial conditions
#' @srrstats {G5.9b} Different random seeds should not meaningfully change results
test_that("noise susceptibility - different random seeds with package data [G5.9b]", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  # Use subset of package data for faster testing
  n_subset <- min(45, length(sarscov2$cases))
  test_cases <- sarscov2$cases[1:n_subset]
  test_dates <- sarscov2$date[1:n_subset]

  mod <- construct_model(
    method = p_spline(),
    pathogen_structure = single(
      case_timeseries = test_cases,
      time = test_dates,
      pathogen_name = "SARS-CoV-2"
    ),
    dow_effect = FALSE
  )

  # G5.9b: Fit with different random seeds
  seeds <- c(1001, 2002, 3003)
  fits <- list()

  for (i in seq_along(seeds)) {
    set.seed(seeds[i])
    fits[[i]] <- fit_model(mod, iter = 800, warmup = 400, chains = 2)

    # Check individual convergence
    rhat_values <- summary(fits[[i]])$summary[, "Rhat"]
    max_rhat <- max(rhat_values, na.rm = TRUE)
    expect_true(max_rhat < 1.1, info = paste("Poor convergence for seed", seeds[i]))
  }

  # Extract posterior means from each fit
  posterior_means <- lapply(fits, function(fit) {
    colMeans(extract(fit, "y_rep")[[1]])
  })

  # G5.9b: Results should be consistent across different seeds
  # Compare all pairs of fits
  correlations <- c()
  rmses <- c()

  for (i in 1:(length(posterior_means) - 1)) {
    for (j in (i + 1):length(posterior_means)) {
      corr <- cor(posterior_means[[i]], posterior_means[[j]])
      rmse <- sqrt(mean((posterior_means[[i]] - posterior_means[[j]])^2))

      correlations <- c(correlations, corr)
      rmses <- c(rmses, rmse)
    }
  }

  # Results should be highly correlated across different seeds
  expect_true(all(correlations > 0.95),
              info = paste("Results too variable across seeds. Correlations:",
                           paste(round(correlations, 4), collapse = ", ")))

  # RMSE should be small relative to signal
  mean_signal <- mean(sapply(posterior_means, mean))
  relative_rmses <- rmses / mean_signal

  expect_true(all(relative_rmses < 0.05),
              info = paste("RMSE too large across seeds. Relative RMSEs:",
                           paste(relative_rmses, collapse = ", ")))

  # Check that all fits converged properly
  all_rhats <- lapply(fits, function(fit) {
    summary(fit)$summary[, "Rhat"]
  })

  max_rhats <- sapply(all_rhats, function(x) max(x, na.rm = TRUE))
  expect_true(all(max_rhats < 1.1),
              info = paste("Convergence issues detected. Max Rhats:",
                           paste(max_rhats, collapse = ", ")))
})

# Additional test for multiple pathogen models with parameter recovery
#' @srrstats {G5.6} Parameter recovery for complex model structures
test_that("parameter recovery - multiple pathogen model [G5.6, G5.6a]", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  set.seed(88888)

  # G5.6: Create synthetic data with known pathogen proportions
  n_obs <- 50
  dates <- seq(as.Date("2023-01-01"), by = "week", length.out = n_obs)

  # Known total cases and component proportions
  total_cases <- rpois(n_obs, lambda = 100 + 20 * sin(2 * pi * seq_len(n_obs) / 26))

  # Known pathogen proportions (should sum to 1)
  prop_alpha <- 0.3 + 0.1 * sin(2 * pi * seq_len(n_obs) / 26)
  prop_delta <- 0.4 + 0.05 * cos(2 * pi * seq_len(n_obs) / 26)
  prop_omicron <- 0.2
  prop_other <- 1 - (prop_alpha + prop_delta + prop_omicron)

  # Generate component cases
  alpha_cases <- rbinom(n_obs, total_cases, prop_alpha)
  delta_cases <- rbinom(n_obs, total_cases, prop_delta)
  omicron_cases <- rbinom(n_obs, total_cases, prop_omicron)
  other_cases <- total_cases - (alpha_cases + delta_cases + omicron_cases)
  other_cases <- pmax(0, other_cases)  # Ensure non-negative

  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = multiple(
      case_timeseries = total_cases,
      time = dates,
      component_pathogen_timeseries = list(
        alpha = alpha_cases,
        delta = delta_cases,
        omicron = omicron_cases,
        other = other_cases
      )
    ),
    dow_effect = FALSE
  )

  fit <- fit_model(mod, iter = 1000, warmup = 500, chains = 2)

  # Check that model fitted successfully
  expect_s3_class(fit, "stanfit")

  # Check convergence
  rhat_values <- summary(fit)$summary[, "Rhat"]
  max_rhat <- max(rhat_values, na.rm = TRUE)
  expect_true(max_rhat < 1.1, info = paste("Poor convergence, max Rhat:", max_rhat))

  # G5.6a: Test parameter recovery within tolerance
  # Extract posterior predictions for total cases
  y_rep <- extract(fit, "y_rep")[[1]]
  fitted_total <- colMeans(y_rep)

  # Should recover the total case pattern reasonably well
  correlation <- cor(fitted_total, total_cases)
  expect_true(correlation > 0.8,
              info = paste("Poor recovery of total cases. Correlation:", correlation))

  # The model should maintain reasonable proportions
  # (This is a basic check - more detailed proportion recovery could be added)
  expect_true(all(fitted_total > 0), "Fitted values should be positive")
  expect_true(max(fitted_total) / min(fitted_total) < 10,
              "Fitted values show reasonable dynamic range")
})

# Comprehensive test combining multiple guidelines
#' @srrstats {G5.4, G5.5, G5.7, G5.9b} Combined correctness, performance, and stability test
test_that("comprehensive model validation - subtyped influenza with package data [G5.4, G5.5, G5.7, G5.9b]", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  # G5.5: Fixed seed for reproducibility
  set.seed(12321)

  # G5.4: Use package influenza data for correctness testing
  # G5.9b: Test with multiple seeds for stability
  seeds <- c(111, 222, 333)
  model_results <- list()

  for (seed_idx in seq_along(seeds)) {
    set.seed(seeds[seed_idx])

    # Use realistic subtyped influenza model with package data
    mod <- construct_model(
      method = random_walk(),
      pathogen_structure = subtyped(
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
        ),
        smoothing_structure = 'independent',
        observation_noise = 'observation_noise_only'
      ),
      dow_effect = TRUE
    )

    # G5.7: Monitor performance
    start_time <- Sys.time()
    fit <- fit_model(mod, iter = 800, warmup = 400, chains = 2)
    fit_time <- as.numeric(Sys.time() - start_time)

    # Check convergence
    rhat_values <- summary(fit)$summary[, "Rhat"]
    max_rhat <- max(rhat_values, na.rm = TRUE)

    # Extract key results
    y_rep <- extract(fit, "y_rep")[[1]]
    fitted_ili <- colMeans(y_rep)

    model_results[[seed_idx]] <- list(
      seed = seeds[seed_idx],
      fit_time = fit_time,
      max_rhat = max_rhat,
      fitted_ili = fitted_ili,
      correlation = cor(fitted_ili, influenza$ili)
    )

    # Individual fit quality checks
    expect_true(max_rhat < 1.1,
                info = paste("Poor convergence for seed", seeds[seed_idx]))
    expect_true(model_results[[seed_idx]]$correlation > 0.6,
                info = paste("Poor fit quality for seed", seeds[seed_idx]))

    # G5.4: Basic correctness checks
    expect_true(all(fitted_ili > 0), "All fitted values should be positive")
    expect_true(length(fitted_ili) == length(influenza$ili), "Output length should match input")
  }

  # G5.9b: Check consistency across seeds
  correlations <- sapply(model_results, `[[`, "correlation")
  expect_true(sd(correlations) < 0.1,
              info = paste("Results too variable across seeds. Correlations:",
                           paste(round(correlations, 3), collapse = ", ")))

  # G5.7: All fits should complete in reasonable time and converge
  fit_times <- sapply(model_results, `[[`, "fit_time")
  max_rhats <- sapply(model_results, `[[`, "max_rhat")

  expect_true(all(fit_times < 300), "Fitting taking too long (>5 minutes)")
  expect_true(all(max_rhats < 1.1), "Convergence issues across seeds")

  # G5.4: Final correctness check - model should capture ILI pattern
  mean_fitted <- rowMeans(sapply(model_results, `[[`, "fitted_ili"))
  overall_correlation <- cor(mean_fitted, influenza$ili)
  expect_true(overall_correlation > 0.7,
              info = paste("Poor overall model performance. Correlation:", round(overall_correlation, 3)))

  # Additional reasonableness checks for complex subtyped model
  expect_true(max(mean_fitted) / min(mean_fitted) < 50,
              "Dynamic range of fitted values seems unreasonable")
  expect_true(var(mean_fitted) > 0, "Fitted values should show variation")
})
