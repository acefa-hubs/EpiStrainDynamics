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
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
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

  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling for all numeric
  #'   inputs
  case_timeseries <- as.numeric(case_timeseries)

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  validate_matching_lengths(
    case_timeseries = case_timeseries,
    time = time
  )

  #' @srrstats {G2.2} restrict to selection of only first provided name
  if (length(pathogen_name) != 1) {
    single_pathogen_name <- pathogen_name[1]
    cli::cli_alert('More than one {pathogen_name} provided, using only
                   {.var {single_pathogen_name')
  }

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
#'   smoothing structure; \code{tau[1]}), `independent` (each pathogen with
#'   independent smoothing structure; \code{tau[number of pathogens]}), or
#'   `correlated` (smoothing structure is correlated among pathogens
#'   \code{Sigma[number of pathogens, number of pathogens]}). Case-insensitive.
#' @param observation_noise either `observation_noise_only` (only includes
#'   observation noise - the same between pathogens) or `pathogen_specific_noise`
#'   (includes noise in individual pathogens as well). Case-insensitive.
#'
#' @returns named list including pathogen_structure, pathogen_names, data,
#'   and model_params of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
#' @export
#'
#' @srrstats {G1.3} clear definitions of smoothing_structure and
#'   observation_noise
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
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

  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling for all numeric
  #'   inputs
  case_timeseries <- as.numeric(case_timeseries)
  component_pathogen_timeseries <- lapply(component_pathogen_timeseries, as.numeric)

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.1a, G2.6} validation for vector inputs
  validate_matching_lengths(
    case_timeseries = case_timeseries,
    time = time
  )
  validate_list_vector(
    component_pathogen_timeseries,
    "component_pathogen_timeseries",
    reference_length = length(case_timeseries)
  )

  #' @srrstats {G2.3, G2.3a, G2.3b} univariate variables must match and be
  #'   case insensitive
  smoothing_structure <- tolower(rlang::arg_match(
    arg = smoothing_structure
  ))
  observation_noise <- tolower(rlang::arg_match(
    arg = observation_noise
  ))
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
#'   smoothing structure; \code{tau[1]}), `independent` (each pathogen with
#'   independent smoothing structure; \code{tau[number of pathogens]}), or
#'   `correlated` (smoothing structure is correlated among pathogens
#'   \code{Sigma[number of pathogens, number of pathogens]}). Case-insensitive.
#' @param observation_noise either `observation_noise_only` (only includes
#'   observation noise - the same between pathogens) or `pathogen_specific_noise`
#'   (includes noise in individual pathogens as well). Case-insensitive.
#'
#' @returns named list including pathogen_structure, pathogen_names, data,
#'   and model_params of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
#' @export
#'
#' @srrstats {G1.3} clear definitions of smoothing_structure and
#'   observation_noise
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
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

  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling for all numeric
  #'   inputs
  case_timeseries <- as.numeric(case_timeseries)
  influenzaA_unsubtyped_timeseries <- as.numeric(influenzaA_unsubtyped_timeseries)
  influenzaA_subtyped_timeseries <- lapply(influenzaA_subtyped_timeseries, as.numeric)
  other_pathogen_timeseries <- lapply(other_pathogen_timeseries, as.numeric)

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.1a, G2.6} validation for vector inputs
  validate_matching_lengths(
    case_timeseries = case_timeseries,
    time = time,
    influenzaA_unsubtyped_timeseries = influenzaA_unsubtyped_timeseries
  )
  validate_list_vector(
    influenzaA_subtyped_timeseries,
    "influenzaA_subtyped_timeseries",
    reference_length = length(case_timeseries)
  )
  validate_list_vector(
    other_pathogen_timeseries,
    "other_pathogen_timeseries",
    reference_length = length(case_timeseries)
  )

  #' @srrstats {G2.3, G2.3a, G2.3b} univariate variables must match and be
  #'   case insensitive
  smoothing_structure <- tolower(rlang::arg_match(
    arg = smoothing_structure
  ))
  observation_noise <- tolower(rlang::arg_match(
    arg = observation_noise
  ))
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
#' @param smoothing_structure either `shared` (all pathogens have the same
#'   smoothing structure; \code{tau[1]}), `independent` (each pathogen with
#'   independent smoothing structure; \code{tau[number of pathogens]}), or
#'   `correlated` (smoothing structure is correlated among pathogens
#'   \code{Sigma[number of pathogens, number of pathogens]}). Case-insensitive.
#'
#' @returns numeric value for covariance structure needed by stan models
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
get_cov_structure <- function(smoothing_structure = c('shared',
                                                      'independent',
                                                      'correlated')) {

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
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
get_noise_structure <- function(observation_noise = c('observation_noise_only',
                                                      'pathogen_specific_noise')) {

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
