# Prevent browser from opening during tests (rstan diagnostics)
options(browser = function(...) {})

# Helper to check package data exists
check_package_data <- function() {
  if (!exists("sarscov2") || !exists("influenza")) {
    skip("Package data (sarscov2, influenza) not available")
  }
}

#' @srrstats {G5.4, G5.5} Correctness tests with fixed test data + fixed seed
test_that("Models produce expected results on fixed test data", {
  skip_on_cran()  # Skip on CRAN due to computational time

  set.seed(98765)

  # Test with sarscov2 single pathogen
  model_single <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      data = sarscov2,
      case_timeseries = 'cases',
      time = 'date'
    )
  )

  fit_single <- fit_model(model_single, n_chain = 3, n_iter = 3000,
                          multi_cores = FALSE, seed = 98765, verbose = FALSE)
  inc_single <- incidence(fit_single, dow = FALSE)

  # Expected values from a known good run (update these after first run)
  # Run once, inspect the output, then hardcode these values
  inc <- inc_single$measure$y
  expect_type(inc, "double")
  expect_length(inc, nrow(sarscov2))
  expect_true(all(inc > 0))
  expect_equal(mean(inc), mean(sarscov2$cases), tolerance = 0.15)

  # Test with multiple pathogens
  model_multi <- construct_model(
    method = random_walk(),
    pathogen_structure = multiple(
      data = sarscov2,
      case_timeseries = 'cases',
      time = 'date',
      component_pathogen_timeseries = c('alpha', 'delta', 'omicron', 'other')
    )
  )

  fit_multi <- fit_model(model_multi, n_chain = 4, n_iter = 3000,
                         multi_cores = FALSE, seed = 98765, verbose = FALSE)
  inc_multi <- incidence(fit_multi, dow = FALSE)
  inc <- inc_multi$measure$y[inc_multi$measure$pathogen == 'Total']

  # Check structural correctness
  expect_type(inc, "double")
  expect_length(inc, nrow(sarscov2))

  # Check total incidence matches
  expect_equal(mean(inc), mean(sarscov2$cases), tolerance = 0.15)
})
