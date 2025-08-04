#' Construct model
#'
#' @param method either random_walk() or p_spline()
#' @param pathogen_structure either single(), multiple(), or subtyped()
#' @param dow_effect logical whether to incorporate a day of week model
#'
#' @returns a list containing the data, the model parameters, and pathogen
#'  names
#' @export
#'
construct_model <- function (method,
                             pathogen_structure,
                             dow_effect = FALSE) {

  model_type <- get_model_type(
    method$method,
    pathogen_structure$pathogen_structure
  )

  time <- seq(1, length(pathogen_structure$data$case_timeseries))

  # Set up day-of-week effects
  if (dow_effect) {
    week_effect <- 7
    DOW <- (time %% 7) + 1
  }
  if (!dow_effect) {
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
    pathogen_names = pathogen_structure$pathogen_names,
    dow_effect = dow_effect
  )

  class(model_input) <- c(model_type, class(model_input))
  return(model_input)
}

#' Single pathogen structure
#'
#' @param case_timeseries total case data
#' @param pathogen_name name of pathogen
#'
#' @returns formatted list with pathogen structure and data
#' @export
#'
single <- function (case_timeseries, pathogen_name) {
  model_inputs <- list(
    pathogen_structure = 'single',
    pathogen_names = pathogen_name,
    data = list(
      case_timeseries = case_timeseries
    )
  )

  class(model_inputs) <- 'EpiStrainDynamics.pathogen_structure'
  model_inputs
}

#' Multiple pathogen structure
#'
#' @param case_timeseries total case data
#' @param component_pathogen_timeseries named list of each pathogen that has
#'   timeseries of case counts
#' @param smoothing_structure either 'shared' (all pathogens have the same;
#'   tau[1]), 'independent' (each pathogen has completely independent smoothing
#'   structure; tau[number of pathogens]), or 'correlated' (smoothing structure is
#'   correlated among pathogens Sigma[number of pathogens, number of pathogens])
#' @param observation_noise either 'observation_noise_only' (only includes
#'   observation noise - the same between pathogens) or pathogen_specific_noise'
#'   (includes noise in individual pathogens as well)
#' @returns named list including pathogen_structure, data, and model_params
#' @export
#'
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

#' Subtyped pathogen structure
#'
#' @param case_timeseries total case data
#' @param influenzaA_unsubtyped_timeseries case timeseries data for unsubtyped
#'   influenzaA
#' @param influenzaA_subtyped_timeseries case timeseries data for subtyped
#'   influenzaA
#' @param other_pathogen_timeseries other pathogen case timeseries
#' @param smoothing_structure either 'shared' (all pathogens have the same;
#'   tau[1]), 'independent' (each pathogen has completely independent smoothing
#'   structure; tau[number of pathogens]), or 'correlated' (smoothing structure is
#'   correlated among pathogens Sigma[number of pathogens, number of pathogens])
#' @param observation_noise either 'observation_noise_only' (only includes
#'   observation noise - the same between pathogens) or pathogen_specific_noise'
#'   (includes noise in individual pathogens as well)
#'
#' @importFrom rlang arg_match
#'
#' @returns named list including pathogen_structure, data, and model_params
#' @export
#'
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

#' Random walk method
#'
#' @returns list with method identified as random walk
#' @export
#'
#' @examples random_walk()
random_walk <- function () {
  out <- list(
    method = 'random-walk'
  )
  return(out)
}

#' P-spline method
#'
#' @param spline_degree polynomial degree of the individual spline segments
#'   used to construct the overall curve
#' @param days_per_knot number of days for each knot
#'
#' @returns list with method and model parameters
#' @export
#'
#' @examples p_spline()
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

#' Extract model type from method and pathogen structure
#'
#' @param method method function, either random_walk() or p_spline()
#' @param pathogens pathogen structure, either single, multiple, or subtyped
#'
#' @returns model type
#'
get_model_type <- function (method, pathogens) {

  is_ps_single <- method == 'p-spline' & pathogens == 'single'
  is_rw_single <- method == 'random-walk' & pathogens == 'single'
  is_ps_multiple <- method == 'p-spline' & pathogens == 'multiple'
  is_rw_multiple <- method == 'random-walk' & pathogens == 'multiple'
  is_ps_subtyped <- method == 'p-spline' & pathogens == 'subtyped'
  is_rw_subtyped <- method == 'random-walk' & pathogens == 'subtyped'

  if (is_ps_single) model_type <- 'ps_single'
  if (is_rw_single) model_type <- 'rw_single'
  if (is_ps_multiple) model_type <- 'ps_multiple'
  if (is_rw_multiple) model_type <- 'rw_multiple'
  if (is_ps_subtyped) model_type <- 'ps_subtyped'
  if (is_rw_subtyped) model_type <- 'rw_subtyped'

  return(model_type)
}

#' Get covariate structure from assigned smoothing structure
#'
#' @param smoothing_structure either 'shared' (all pathogens have the same;
#'   tau[1]), 'independent' (each pathogen has completely independent smoothing
#'   structure; tau[number of pathogens]), or 'correlated' (smoothing structure is
#'   correlated among pathogens Sigma[number of pathogens, number of pathogens])
#'
#' @returns numeric value for covariance structure needed by stan models
#'
get_cov_structure <- function (smoothing_structure) {
  if (smoothing_structure == 'single') cov_structure <- 0
  if (smoothing_structure == 'independent') cov_structure <- 1
  if (smoothing_structure == 'correlated') cov_structure <- 2
  return(cov_structure)
}

#' Get noise structure from specified observation noise
#'
#' @param observation_noise either 'observation_noise_only' (only includes
#'   observation noise - the same between pathogens) or pathogen_specific_noise'
#'   (includes noise in individual pathogens as well)
#'
#' @returns numeric value for noise structure needed by stan models
#'
get_noise_structure <- function (observation_noise) {
  if (observation_noise == 'observation_noise_only') noise_structure <- 0
  if (observation_noise == 'pathogen_specific_noise') noise_structure <- 1
  return(noise_structure)
}
