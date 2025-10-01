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

  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  if (any(is.na(case_timeseries)) || any(is.na(time))) {
    cli::cli_abort("{case_timeseries} and {time} cannot contain NA values.
                   `EpiStrainDynamics` expects data in regular timesteps with
                   no gaps. Please review your data.")
  }

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
#'
#' @returns named list including pathogen_structure, pathogen_names, and data
#'   of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
#' @export
#'
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
#'   )
#' )
#'
multiple <- function (case_timeseries,
                      time,
                      component_pathogen_timeseries) {

  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  na_in_components <- do.call(
    any, lapply(component_pathogen_timeseries, function(x) is.na(x))
  )
  if (any(is.na(case_timeseries)) || any(is.na(time)) || na_in_components) {
    cli::cli_abort("Input data cannot contain NA values. `EpiStrainDynamics`
    expects data in regular timesteps with no gaps. Please review your data.")
  }

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
#'
#' @returns named list including pathogen_structure, pathogen_names, and data
#'   of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
#' @export
#'
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
#'   )
#' )
#'
subtyped <- function (case_timeseries,
                      time,
                      influenzaA_unsubtyped_timeseries,
                      influenzaA_subtyped_timeseries,
                      other_pathogen_timeseries) {

  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  na_unsubtyped <- do.call(
    any, lapply(influenzaA_unsubtyped_timeseries, function(x) is.na(x))
  )
  na_subtyped <- do.call(
    any, lapply(influenzaA_subtyped_timeseries, function(x) is.na(x))
  )
  na_other <- do.call(
    any, lapply(other_pathogen_timeseries, function(x) is.na(x))
  )
  if (any(is.na(case_timeseries)) || any(is.na(time)) || na_unsubtyped
      || na_subtyped || na_other) {
    cli::cli_abort("Input data cannot contain NA values. `EpiStrainDynamics`
    expects data in regular timesteps with no gaps. Please review your data.")
  }

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
    )
  )

  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}
