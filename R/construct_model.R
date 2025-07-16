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

#' Single pathogen structure
#'
#' @param case_timeseries total case data
#'
#' @returns formatted list with pathogen structure and data
#' @export
#'
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

