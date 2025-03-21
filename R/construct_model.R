# this approach leaves less room for a user to put in the wrong data or select options that are irrelevant.
# makes the choies they should make more explicit rather than needing to follow branching paths of
# documentation perfectly and us relying on lots of stop>errors or messages.

# also separates the decisions needed for data prep and decisions needed to model parameterisation.
construct_model <- function (method,
                             pathogen_structure,
                             dow_data = NULL) {

  model_type <- get_model_type(
    method$method,
    pathogen_structure$pathogen_structure
  )

  # Set up day-of-week effects
  if (!is.null(dow_data)) {
    week_effect <- 7
    DOW <- (dow_data %% 7) + 1
  }
  if (is.null(dow_data)) {
    time <- seq(1, length(pathogen_structure$data$total_cases))
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
      list(knots = method$knots)
    )
  }

  model_input <- list(
    data = pathogen_structure$data,
    model_params = model_params,
    pathogen_names = pathogen_structure$pathogen_names
  )

  class(model_input) <- c(model_type, class(model_input))
  return(model_input)
}




single <- function (total_case_data) {
  model_inputs <- list(
    pathogen_structure = 'single',
    data = list(
      total_cases = total_case_data
    )
  )

  class(model_inputs) <- 'EpiStrainDynamics.pathogen_structure'
  model_inputs
}

multiple <- function (total_case_data,
                      component_pathogen_data,
                      smoothing_structure = c(
                        'single',
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

  pathogen_names <- names(component_pathogen_data)

  component_pathogens <- t(matrix(
    data = do.call(c, component_pathogen_data),
    ncol = length(component_pathogen_data)
  ))

  model_inputs <- list(
    pathogen_structure = 'multiple',
    pathogen_names = pathogen_names,
    data = list(
      total_cases = total_case_data,
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

subtyped <- function (total_case_data,
                      influenzaA_data,
                      influenzaA_subtype_data,
                      other_component_pathogen_data,
                      smoothing_structure = c(
                        'single',
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

  # add check that influenzaA_subtype_data and other_component_pathogen_data are named lists
  pathogen_names <- c(
    names(influenzaA_subtype_data),
    names(other_component_pathogen_data)
  )

  component_pathogens_list <- append(
    list(influenzaA = influenzaA_data),
    other_component_pathogen_data
  )

  component_pathogens <- t(matrix(
    data = do.call(c, component_pathogens_list),
    ncol = length(component_pathogens_list)
  ))
  influenzaA_subtypes <- t(matrix(
    data = do.call(c, influenzaA_subtype_data),
    ncol = length(influenzaA_subtype_data)
  ))

  model_inputs <- list(
    pathogen_structure = 'subtyped',
    pathogen_names = pathogen_names,
    data = list(
      total_cases = total_case_data,
      component_pathogens = component_pathogens,
      influenzaA_subtypes = influenzaA_subtypes
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

p_spline <- function (time_data,
                      spline_degree = 3,
                      days_per_knot = 3) {

  # Calculate the locations of equally spaced knots
  knots <- get_knots(time_data,
                     days_per_knot = days_per_knot,
                     spline_degree = spline_degree)

  out <- list(
    method = 'p-spline',
    knots = knots
  )
  return(out)
}

