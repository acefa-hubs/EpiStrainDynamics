# this approach leaves less room for a user to put in the wrong data or select options that are irrelevant.
# makes the choies they should make more explicit rather than needing to follow branching paths of
# documentation perfectly and us relying on lots of stop>errors or messages.

# also separates the decisions needed for data prep and decisions needed to model parameterisation.
construct_model <- function (method,
                             pathogen_structure,
                             dow_effect = NULL) {

  model_type <- get_model_type(
    method$method,
    pathogen_structure$pathogen_structure
  )

  time <- seq(1, length(pathogen_structure$data$case_timeseries))

  # Set up day-of-week effects
  if (!is.null(dow_effect)) {
    week_effect <- 7
    DOW <- (time %% 7) + 1
  }
  if (is.null(dow_effect)) {
    week_effect <- 1
    DOW <- (time %% 1) + 1
  }

  model_params <- append(
    list(
      week_effect = week_effect,
      DOW = DOW
    ),
    pathogen_structure$model_params)

  if (method$method == 'p-spline') {
    model_params <- append(
      model_params,
      method$model_params
    )
  }

  data <- append(
    list(time = time),
    pathogen_structure$data
  )

  model_input <- list(
    data = data,
    model_params = model_params,
    pathogen_names = pathogen_structure$pathogen_names
  )

  class(model_input) <- c(model_type, class(model_input))
  return(model_input)
}

single <- function (case_timeseries) {
  model_inputs <- list(
    pathogen_structure = 'single',
    data = list(
      case_timeseries = case_timeseries
    )
  )

  class(model_inputs) <- 'EpiStrainDynamics.pathogen_structure'
  model_inputs
}

multiple <- function (case_timeseries,
                      component_pathogen_timeseries,
                      smoothing_structure = c(
                        'shared',
                        'independent',
                        'correlated'),
                      observation_noise = c(
                        'observation_noise_only',
                        'pathogen_specific_noise')) {

  smoothing_structure <- rlang::arg_match(
    arg = smoothing_structure
  )
  observation_noise <- rlang::arg_match(
    arg = observation_noise
  )

  cov_structure <- get_cov_structure(smoothing_structure)
  noise_structure <- get_noise_structure(observation_noise)

  pathogen_names <- names(component_pathogen_timeseries)

  component_pathogens <- t(matrix(
    data = do.call(c, component_pathogen_timeseries),
    ncol = length(component_pathogen_timeseries)
  ))

  model_inputs <- list(
    pathogen_structure = 'multiple',
    pathogen_names = pathogen_names,
    data = list(
      case_timeseries = case_timeseries,
      component_pathogens = component_pathogens
    ),
    model_params = list(
      cov_structure = cov_structure,
      noise_structure = noise_structure
    )
  )

  class(model_inputs) <- 'EpiStrainDynamics.pathogen_structure'
  model_inputs
}

subtyped <- function (case_timeseries,
                      influenzaA_unsubtyped_timeseries,
                      influenzaA_subtyped_timeseries,
                      other_pathogen_timeseries,
                      smoothing_structure = c(
                        'shared',
                        'independent',
                        'correlated'),
                      observation_noise = c(
                        'observation_noise_only',
                        'pathogen_specific_noise')) {

  smoothing_structure <- rlang::arg_match(
    arg = smoothing_structure
  )
  observation_noise <- rlang::arg_match(
    arg = observation_noise
  )

  cov_structure <- get_cov_structure(smoothing_structure)
  noise_structure <- get_noise_structure(observation_noise)

  # add check that influenzaA_subtyped_timeseries and other_pathogen_timeseries are named lists
  pathogen_names <- c(
    names(influenzaA_subtyped_timeseries),
    names(other_pathogen_timeseries)
  )

  component_pathogens_list <- append(
    list(influenzaA = influenzaA_unsubtyped_timeseries),
    other_pathogen_timeseries
  )

  component_pathogens <- t(matrix(
    data = do.call(c, component_pathogens_list),
    ncol = length(component_pathogens_list)
  ))
  influenzaA_subtyped <- t(matrix(
    data = do.call(c, influenzaA_subtyped_timeseries),
    ncol = length(influenzaA_subtyped_timeseries)
  ))

  model_inputs <- list(
    pathogen_structure = 'subtyped',
    pathogen_names = pathogen_names,
    data = list(
      case_timeseries = case_timeseries,
      component_pathogens = component_pathogens,
      influenzaA_subtyped = influenzaA_subtyped
    ),
    model_params = list(
      cov_structure = cov_structure,
      noise_structure = noise_structure
    )
  )

  class(model_inputs) <- 'EpiStrainDynamics.pathogen_structure'
  model_inputs
}

random_walk <- function () {
  out <- list(
    method = 'random-walk'
  )
  return(out)
}

p_spline <- function (spline_degree = 3,
                      days_per_knot = 3) {

  out <- list(
    method = 'p-spline',
    model_params = list(
      spline_degree = spline_degree,
      days_per_knot = days_per_knot
    )
  )
  return(out)
}

