# Files that use these: test-fit-model.R, test-diagnose-model.R,
# test-metrics.R, test-plot.R

# Only create fitted models if package data is available
if (exists("sarscov2") && exists("influenza")) {

  message("Setting up cached fitted models for test session...")

  # Suppress Stan output during setup
  options(browser = function(...) {})

  # Small subsets for speed
  sarscov2_subset <- sarscov2[1:50, ]
  influenza_subset <- influenza[1:40, ]

  # Create and cache fitted models in the global test environment
  # These will be available to ALL test files

  fit_rw_single <<- suppressWarnings(
    fit_model(
      construct_model(
        method = random_walk(),
        pathogen_structure = single(
          data = sarscov2_subset,
          case_timeseries = 'cases',
          time = 'date'
        )
      ),
      n_iter = 500,
      n_chain = 1,
      verbose = FALSE,
      seed = 123
    )
  )

  fit_ps_single <<- suppressWarnings(
    fit_model(
      construct_model(
        method = p_spline(),
        pathogen_structure = single(
          data = sarscov2_subset,
          case_timeseries = 'cases',
          time = 'date'
        )
      ),
      n_iter = 500,
      n_chain = 1,
      verbose = FALSE,
      seed = 123
    )
  )

  fit_rw_multi <<- suppressWarnings(
    fit_model(
      construct_model(
        method = random_walk(),
        pathogen_structure = multiple(
          data = sarscov2_subset,
          case_timeseries = 'cases',
          time = 'date',
          component_pathogen_timeseries = c('alpha', 'delta',
                                            'omicron', 'other')
        )
      ),
      n_iter = 500,
      n_chain = 1,
      verbose = FALSE,
      seed = 123
    )
  )

  fit_ps_multi <<- suppressWarnings(
    fit_model(
      construct_model(
        method = p_spline(),
        pathogen_structure = multiple(
          data = sarscov2_subset,
          case_timeseries = 'cases',
          time = 'date',
          component_pathogen_timeseries = c('alpha', 'delta',
                                            'omicron', 'other')
        )
      ),
      n_iter = 500,
      n_chain = 1,
      verbose = FALSE,
      seed = 123
    )
  )

  fit_rw_subtyped <<- suppressWarnings(
    fit_model(
      construct_model(
        method = random_walk(),
        pathogen_structure = subtyped(
          data = influenza_subset,
          case_timeseries = 'ili',
          time = 'week',
          influenzaA_unsubtyped_timeseries = 'inf_A',
          influenzaA_subtyped_timeseries = c('inf_H3N2', 'inf_H1N1'),
          other_pathogen_timeseries = c('inf_B', 'other')
        )
      ),
      n_iter = 500,
      n_chain = 1,
      verbose = FALSE,
      seed = 123
    )
  )

  fit_ps_subtyped <<- suppressWarnings(
    fit_model(
      construct_model(
        method = p_spline(),
        pathogen_structure = subtyped(
          data = influenza_subset,
          case_timeseries = 'ili',
          time = 'week',
          influenzaA_unsubtyped_timeseries = 'inf_A',
          influenzaA_subtyped_timeseries = c('inf_H3N2', 'inf_H1N1'),
          other_pathogen_timeseries = c('inf_B', 'other')
        )
      ),
      n_iter = 500,
      n_chain = 1,
      verbose = FALSE,
      seed = 123
    )
  )

  fit_rw_single_dow <<- suppressWarnings(
    fit_model(
      construct_model(
        method = random_walk(),
        pathogen_structure = single(
          data = sarscov2_subset,
          case_timeseries = 'cases',
          time = 'date'
        ),
        dow_effect = TRUE
      ),
      n_iter = 500,
      n_chain = 1,
      verbose = FALSE,
      seed = 123
    )
  )

  # Helper function for generation interval
  gi_simple <<- function(x) {
    ifelse(x == 0, 0, 4 * x * exp(-2 * x))
  }

  inc_single <<- incidence(fit_rw_single, dow = FALSE)
  inc_multi <<- incidence(fit_rw_multi, dow = FALSE)
  gr_single <<- growth_rate(fit_rw_single)
  gr_multi <<- growth_rate(fit_rw_multi)
  rt_single <<- Rt(fit_rw_single, tau_max = 7, gi_dist = gi_simple)
  rt_multi <<- Rt(fit_rw_multi, tau_max = 7, gi_dist = gi_simple)
  prop <<- proportion(fit_rw_multi)

  message("Cached fitted models ready: fit_rw_single, fit_ps_single, fit_rw_multi, fit_ps_multi, fit_rw_subtyped, fit_ps_subtyped, fit_rw_single_dow. Cached metrics outputs ready: inc_single, inc_multi, gr_single, gr_multi, rt_single, rt_multi, prop")

} else {
  message("Package data not available - skipping fitted model setup")
}
