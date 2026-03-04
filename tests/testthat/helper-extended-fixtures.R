# tests/testthat/helper-extended-fixtures.R
#
# Generates fitted model fixtures for extended tests.
# This file is sourced automatically by testthat before tests run.
# The fixtures are only generated when extended tests are enabled via:
#   EPISTRAINDYNAMICS_EXTENDED_TESTS="true"
#
# Uses real package data (sarscov2, influenza) with more thorough MCMC settings
# (higher iterations, more chains) to validate models work correctly on actual
# surveillance data patterns.

extended_tests_enabled <- function() {
  identical(Sys.getenv("EPISTRAINDYNAMICS_EXTENDED_TESTS"), "true")
}

skip_extended <- function() {
  if (!extended_tests_enabled()) {
    skip("Extended tests disabled. Set EPISTRAINDYNAMICS_EXTENDED_TESTS=true to run.")
  }
}

# This environment will hold the fixtures so tests can access them
# without re-fitting models. Populated only when extended tests are on.
extended_fixtures <- new.env(parent = emptyenv())

if (extended_tests_enabled()) {

  message("EpiStrainDynamics: generating extended test fixtures (this may take a few minutes)...")

  # Check that package data is available
  if (!exists("sarscov2") || !exists("influenza")) {
    stop("Package data (sarscov2, influenza) not available for extended tests")
  }

  # --------------------------------------------------------------------------
  # Store package data in fixtures for test access
  # --------------------------------------------------------------------------
  extended_fixtures$sarscov2 <- sarscov2
  extended_fixtures$influenza <- influenza

  # --------------------------------------------------------------------------
  # Fit models with thorough MCMC settings (more iterations, full chains)
  # --------------------------------------------------------------------------

  # -- Random walk, single pathogen (sarscov2 total cases) --
  message("  Fitting: Random walk, single pathogen...")
  extended_fixtures$fit_rw_single <- fit_model(
    construct_model(
      method = random_walk(),
      pathogen_structure = single(
        data = sarscov2,
        case_timeseries = "cases",
        time = "date"
      )
    ),
    n_iter = 2000,
    n_chain = 4,
    verbose = FALSE
  )

  # -- Random walk, single pathogen with day-of-week effects --
  # message("  Fitting: Random walk, single pathogen, day-of-week...")
  # extended_fixtures$fit_rw_single_dow <- fit_model(
  #   construct_model(
  #     method = random_walk(),
  #     pathogen_structure = single(
  #       data = sarscov2,
  #       case_timeseries = "cases",
  #       time = "date"
  #     ),
  #     dow_effect = TRUE
  #   ),
  #   n_iter = 2000,
  #   n_chain = 4,
  #   verbose = FALSE
  # )

  # -- P-spline, single pathogen --
  message("  Fitting: P-spline, single pathogen...")
  extended_fixtures$fit_ps_single <- fit_model(
    construct_model(
      method = p_spline(),
      pathogen_structure = single(
        data = sarscov2,
        case_timeseries = "cases",
        time = "date"
      )
    ),
    n_iter = 2000,
    n_chain = 4,
    verbose = FALSE
  )

  # -- Random walk, multiple pathogens (sarscov2 variants) --
  message("  Fitting: Random walk, multiple pathogens...")
  extended_fixtures$fit_rw_multi <- fit_model(
    construct_model(
      method = random_walk(),
      pathogen_structure = multiple(
        data = sarscov2,
        case_timeseries = "cases",
        time = "date",
        component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
      )
    ),
    n_iter = 2000,
    n_chain = 4,
    verbose = FALSE
  )

  # -- P-spline, multiple pathogens --
  message("  Fitting: P-spline, multiple pathogens...")
  extended_fixtures$fit_ps_multi <- fit_model(
    construct_model(
      method = p_spline(),
      pathogen_structure = multiple(
        data = sarscov2,
        case_timeseries = "cases",
        time = "date",
        component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
      )
    ),
    n_iter = 2000,
    n_chain = 4,
    verbose = FALSE
  )

  # -- Random walk, subtyped pathogens (influenza) --
  message("  Fitting: Random walk, subtyped pathogens...")
  extended_fixtures$fit_rw_subtyped <- fit_model(
    construct_model(
      method = random_walk(),
      pathogen_structure = subtyped(
        data = influenza,
        case_timeseries = "ili",
        time = "week",
        influenzaA_unsubtyped_timeseries = "inf_A",
        influenzaA_subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
        other_pathogen_timeseries = c("inf_B", "other")
      )
    ),
    n_iter = 2000,
    n_chain = 4,
    verbose = FALSE
  )

  # -- P-spline, subtyped pathogens --
  message("  Fitting: P-spline, subtyped pathogens...")
  extended_fixtures$fit_ps_subtyped <- fit_model(
    construct_model(
      method = p_spline(),
      pathogen_structure = subtyped(
        data = influenza,
        case_timeseries = "ili",
        time = "week",
        influenzaA_unsubtyped_timeseries = "inf_A",
        influenzaA_subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
        other_pathogen_timeseries = c("inf_B", "other")
      )
    ),
    n_iter = 2000,
    n_chain = 4,
    verbose = FALSE
  )

  message("EpiStrainDynamics: extended test fixtures ready.")
}
