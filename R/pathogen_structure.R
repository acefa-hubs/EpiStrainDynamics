#' Single pathogen structure
#'
#' @param case_timeseries vector of total case data
#' @param time vector of labelled time data
#' @param pathogen_name optionally provided name of pathogen. This value is
#'   preserved exactly as provided (case-sensitive). Default name is `default`.
#'
#' @returns formatted list with pathogen structure and data of class
#'   `EpiStrainDynamics.pathogen_structure`.
#' @family pathogen_structure
#' @export
#'
#' @examples
#' single(
#'   case_timeseries = sarscov2$cases,
#'   time = sarscov2$date,
#'   pathogen_name = 'SARS-COV-2' # preserved as provided
#' )
#'
#' single(
#'   case_timeseries = sarscov2$cases,
#'   time = sarscov2$date,
#'   pathogen_name = 'sarscov2' # preserved as provided
#' )
#'
#' single(
#'   case_timeseries = sarscov2$cases,
#'   time = sarscov2$date
#' )
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

  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}

#' Multiple pathogen structure
#'
#' @param case_timeseries vector of total case data
#' @param time vector of labelled time data
#' @param component_pathogen_timeseries named list of each pathogen that has
#'   timeseries of case counts
#' @param smoothing_structure either `shared` (all pathogens have the same
#'   smoothing structure; tau[1]), `independent` (each pathogen has completely
#'   independent smoothing structure; tau[number of pathogens]), or `correlated`
#'   (smoothing structure is correlated among pathogens Sigma[number of
#'   pathogens, number of pathogens]). Case-insensitive.
#' @param observation_noise either `observation_noise_only` (only includes
#'   observation noise - the same between pathogens) or `pathogen_specific_noise`
#'   (includes noise in individual pathogens as well). Case-insensitive.
#'
#' @returns named list including pathogen_structure, pathogen_names, data,
#'   and model_params of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
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
  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
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
#' @param smoothing_structure either `shared` (all pathogens have the same
#'   smoothing structure; tau[1]), `independent` (each pathogen has completely
#'   independent smoothing structure; tau[number of pathogens]), or `correlated`
#'   (smoothing structure is correlated among pathogens Sigma[number of
#'   pathogens, number of pathogens]). Case-insensitive.
#' @param observation_noise either `observation_noise_only` (only includes
#'   observation noise - the same between pathogens) or `pathogen_specific_noise`
#'   (includes noise in individual pathogens as well). Case-insensitive.
#'
#' @returns named list including pathogen_structure, pathogen_names, data,
#'   and model_params of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
#' @export
#'
#' @examples
#' # Both integer and numeric inputs work consistently
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
#' # NOTE: smoothing_structure and observation_noise are case insensitive, so
#' # these will work
#' # subtyped(
#' #    ...
#' #    smoothing_structure = 'Independent',
#' #    observation_noise = 'OBSERVATION_NOISE_ONLY'
#' # )
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

  # Ensure consistent numeric handling for all numeric inputs
  case_timeseries <- as.numeric(case_timeseries)
  influenzaA_unsubtyped_timeseries <- as.numeric(influenzaA_unsubtyped_timeseries)

  # Handle list inputs consistently
  influenzaA_subtyped_timeseries <- lapply(influenzaA_subtyped_timeseries, as.numeric)
  other_pathogen_timeseries <- lapply(other_pathogen_timeseries, as.numeric)

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

  # Ensure consistent matrix creation
  component_pathogens <- t(matrix(
    data = as.numeric(do.call(c, component_pathogens_list)),
    ncol = length(component_pathogens_list)
  ))

  influenzaA_subtyped <- t(matrix(
    data = as.numeric(do.call(c, influenzaA_subtyped_timeseries)),
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

  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}

#' Get covariance structure from assigned smoothing structure
#'
#' @param smoothing_structure either 'shared' (all pathogens have the same;
#'   tau[1]), 'independent' (each pathogen has completely independent smoothing
#'   structure; tau[number of pathogens]), or 'correlated' (smoothing structure is
#'   correlated among pathogens Sigma[number of pathogens, number of pathogens]).
#'   Case-insensitive.
#'
#' @returns numeric value for covariance structure needed by stan models
#'
get_cov_structure <- function(smoothing_structure = c('shared',
                                                      'independent',
                                                      'correlated')) {

  # Handle case insensitivity
  smoothing_structure <- match_arg_case_insensitive(
    arg = smoothing_structure,
    choices = c('shared', 'independent', 'correlated'),
    arg_name = "smoothing_structure"
  )

  # Convert to numeric codes
  cov_structure <- switch(smoothing_structure,
                          'shared' = 0,
                          'independent' = 1,
                          'correlated' = 2,                            {
                            cli::cli_abort("Invalid option provided: '{smoothing_structure}'.
                                             Please choose 'shared', 'independent', or 'correlated'.")
                          }
  )

  return(cov_structure)
}

#' Get noise structure from specified observation noise
#'
#' @param observation_noise either 'observation_noise_only' (only includes
#'   observation noise - the same between pathogens) or 'pathogen_specific_noise'
#'   (includes noise in individual pathogens as well). Case-insensitive.
#'
#' @returns numeric value for noise structure needed by stan models
#'
get_noise_structure <- function(observation_noise = c('observation_noise_only',
                                                      'pathogen_specific_noise')) {

  # Handle case insensitivity
  observation_noise <- match_arg_case_insensitive(
    arg = observation_noise,
    choices = c('observation_noise_only', 'pathogen_specific_noise'),
    arg_name = "observation_noise"
  )

  # Convert to numeric codes
  noise_structure <- switch(observation_noise,
                            'observation_noise_only' = 0,
                            'pathogen_specific_noise' = 1,
                            {
                              cli::cli_abort("Invalid option provided: '{observation_noise}'.
                                             Please choose 'observation_noise_only' or 'pathogen_specific_noise'.")
                            }
  )

  return(noise_structure)
}
