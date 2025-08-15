#' Single pathogen structure
#'
#' @param case_timeseries vector of total case data
#' @param time vector of labelled time data
#' @param pathogen_name name of pathogen
#'
#' @returns formatted list with pathogen structure and data
#' @export
#'
#' @examples
#' single(
#'   case_timeseries = sarscov2$cases,
#'   time = sarscov2$date,
#'   pathogen_name = 'SARS-COV-2'
#' )
#'
single <- function (case_timeseries,
                    time,
                    pathogen_name = 'default') {

  # Input validation
  validate_matching_lengths(
    case_timeseries = case_timeseries,
    time = time
  )

  model_inputs <- list(
    pathogen_structure = 'single',
    pathogen_names = pathogen_name,
    data = list(
      case_timeseries = case_timeseries,
      time = time
    )
  )

  class(model_inputs) <- 'EpiStrainDynamics.pathogen_structure'
  model_inputs
}

#' Multiple pathogen structure
#'
#' @param case_timeseries vector of total case data
#' @param time vector of labelled time data
#' @param component_pathogen_timeseries named list of each pathogen that has
#'   timeseries of case counts
#' @param smoothing_structure either 'shared' (all pathogens have the same
#'   smoothing structure; tau[1]), 'independent' (each pathogen has completely
#'   independent smoothing structure; tau[number of pathogens]), or 'correlated'
#'   (smoothing structure is correlated among pathogens Sigma[number of
#'   pathogens, number of pathogens])
#' @param observation_noise either 'observation_noise_only' (only includes
#'   observation noise - the same between pathogens) or pathogen_specific_noise'
#'   (includes noise in individual pathogens as well)
#' @returns named list including pathogen_structure, pathogen_names, data,
#'   and model_params
#' @export
#'
#' @examples
#' multiple(
#'   case_timeseries = sarscov2$cases,
#'   time = sarscov2$date,
#'   component_pathogen_timeseries = list(
#'     alpha = sarscov2$alpha,
#'     delta = sarscov2$delta,
#'     omicron = sarscov2$omicron,
#'     other = sarscov2$other
#'   ),
#'   smoothing_structure = 'independent',
#'   observation_noise = 'observation_noise_only'
#' )
#'
multiple <- function (case_timeseries,
                      time,
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

  # Input validation
  validate_matching_lengths(
    case_timeseries = case_timeseries,
    time = time
  )

  validate_list_vector_lengths(
    component_pathogen_timeseries,
    "component_pathogen_timeseries",
    reference_length = length(case_timeseries)
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
      time = time,
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
#' @param case_timeseries vector of total case data
#' @param time vector of labelled time data
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
#' @returns named list including pathogen_structure, pathogen_names, data,
#'   and model_params
#' @export
#'
#' @examples
#' subtyped(
#'   case_timeseries = influenza$ili,
#'   time = influenza$week,
#'   influenzaA_unsubtyped_timeseries = influenza$inf_A,
#'   influenzaA_subtyped_timeseries = list(
#'     influenzaA.H3N2 = influenza$inf_H3N2,
#'     influenzaA.H1N1 = influenza$inf_H1N1
#'   ),
#'   other_pathogen_timeseries = list(
#'     influenzaB = influenza$inf_B,
#'     other = influenza$num_spec - influenza$inf_all
#'   ),
#'   smoothing_structure = 'independent',
#'   observation_noise = 'observation_noise_only'
#' )
#'
subtyped <- function (case_timeseries,
                      time,
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

  # Input validation
  validate_matching_lengths(
    case_timeseries = case_timeseries,
    time = time,
    influenzaA_unsubtyped_timeseries = influenzaA_unsubtyped_timeseries
  )

  validate_list_vector_lengths(
    influenzaA_subtyped_timeseries,
    "influenzaA_subtyped_timeseries",
    reference_length = length(case_timeseries)
  )

  validate_list_vector_lengths(
    other_pathogen_timeseries,
    "other_pathogen_timeseries",
    reference_length = length(case_timeseries)
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
      time = time,
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

#' Get covariate structure from assigned smoothing structure
#'
#' @param smoothing_structure either 'shared' (all pathogens have the same;
#'   tau[1]), 'independent' (each pathogen has completely independent smoothing
#'   structure; tau[number of pathogens]), or 'correlated' (smoothing structure is
#'   correlated among pathogens Sigma[number of pathogens, number of pathogens])
#'
#' @returns numeric value for covariance structure needed by stan models
#'
#' @examples
#' get_cov_structure('single')
#'
get_cov_structure <- function (smoothing_structure = c(
                                 'shared',
                                 'independent',
                                 'correlated')) {

  if (smoothing_structure == 'shared') cov_structure <- 0
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
#' @examples
#' get_noise_structure('observation_noise_only')
#'
get_noise_structure <- function (observation_noise = c(
                                   'observation_noise_only',
                                   'pathogen_specific_noise')) {

  if (observation_noise == 'observation_noise_only') noise_structure <- 0
  if (observation_noise == 'pathogen_specific_noise') noise_structure <- 1

  return(noise_structure)
}
