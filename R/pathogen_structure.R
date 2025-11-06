#' Single pathogen structure
#'
#' @param data dataframe containing columns with all relevant data
#' @param case_timeseries name of column with total case data, this name will be
#'   used as the name of the pathogen
#' @param time name of column with time data. flexible format - can be date,
#'   index, or others, accepted as `index` identifiers in the `tsibble` time
#'   format.
#'
#' @returns formatted list with pathogen structure and data of class
#'   `EpiStrainDynamics.pathogen_structure`.
#' @family pathogen_structure
#' @export
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#' @srrstats {G2.7, G2.8} data argument accepts standard tabular forms, which
#'   are then post-processed internally to ensure subsequent functions receive
#'   standard inputs
#'
#' @examples
#' single(
#'   data = sarscov2,
#'   case_timeseries = 'cases',
#'   time = 'date'
#' )
#'
single <- function (data,
                    case_timeseries,
                    time) {

  #' @srrstats {G5.8c, G5.8d} edge cases produce errors
  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.1a} validation for vector inputs
  #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling for all numeric
  #'   inputs
  check_column_exists(data, time, "time")
  check_column_exists(data, case_timeseries, "case_timeseries")
  check_column_numeric(data, case_timeseries, "case_timeseries")

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  all_cols <- c(time, case_timeseries)
  validated_tsbl <- create_validated_tsibble(data, all_cols, time)

  #' @srrstats {G5.3} data objects are returns with no missing values. this
  #'   is validated through conversion to tsibble in `create_validated_tsibble`
  #' @srrstats {G2.10} extraction of single columns happens after conversion
  #'   to consistent tabular form
  model_inputs <- list(
    pathogen_structure = 'single',
    pathogen_names = case_timeseries,
    data = list(
      case_timeseries = validated_tsbl[[case_timeseries]],
      time = validated_tsbl[[time]]
    )
  )

  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}

#' Multiple pathogen structure
#'
#' @param data dataframe containing columns with all relevant data
#' @param case_timeseries name of column with total case data, this name will be
#'   used as the name of the pathogen
#' @param time name of column with time data. flexible format - can be date,
#'   index, or others, accepted as `index` identifiers in the `tsibble` time
#'   format.
#' @param component_pathogen_timeseries vector of column names with additional
#'   pathogen case count timeseries to model
#'
#' @returns named list including pathogen_structure, pathogen_names, and data
#'   of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
#' @export
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#' @srrstats {G2.7, G2.8} data argument accepts standard tabular forms, which
#'   are then post-processed internally to ensure subsequent functions receive
#'   standard inputs
#'
#' @examples
#' multiple(
#'   data = sarscov2,
#'   case_timeseries = 'cases',
#'   time = 'date',
#'   component_pathogen_timeseries = c('alpha', 'delta', 'omicron', 'other')
#' )
#'
multiple <- function (data,
                      case_timeseries,
                      time,
                      component_pathogen_timeseries) {

  #' @srrstats {G5.8c, G5.8d} edge cases produce errors
  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.1a} validation for vector inputs
  #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling for all numeric
  #'   inputs
  check_column_exists(data, case_timeseries, "case_timeseries")
  check_column_numeric(data, case_timeseries, "case_timeseries")
  check_column_exists(data, time, "time")

  for (col in component_pathogen_timeseries) {
    check_column_exists(data, col, "component_pathogen_timeseries")
    check_column_numeric(data, col, "component_pathogen_timeseries")
  }

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  all_cols <- c(time, case_timeseries, component_pathogen_timeseries)
  validated_tsbl <- create_validated_tsibble(data, all_cols, time)

  # Create matrix from validated tsibble
  component_pathogens <- t(as.matrix(validated_tsbl[, component_pathogen_timeseries]))

  #' @srrstats {G5.3} data objects are returns with no missing values. this
  #'   is validated through conversion to tsibble in `create_validated_tsibble`
  #' @srrstats {G2.10} extraction of single columns happens after conversion
  #'   to consistent tabular form
  model_inputs <- list(
    pathogen_structure = 'multiple',
    pathogen_names = component_pathogen_timeseries,
    data = list(
      case_timeseries = validated_tsbl[[case_timeseries]],
      time = validated_tsbl[[time]],
      component_pathogens = component_pathogens
    )
  )
  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}

#' Subtyped pathogen structure
#'
#' @param data dataframe containing columns with all relevant data
#' @param case_timeseries name of column with total case data, this name will be
#'   used as the name of the pathogen
#' @param time name of column with time data. flexible format - can be date,
#'   index, or others, accepted as `index` identifiers in the `tsibble` time
#'   format.
#' @param influenzaA_unsubtyped_timeseries vector of column names with additional
#'   unsubtyped influenzaA case count timeseries
#' @param influenzaA_subtyped_timeseries vector of column names with additional
#'   subtyped influenzaA case count timeseries
#' @param other_pathogen_timeseries vector of column names with additional
#'   pathogen case count timeseries to model
#'
#' @returns named list including pathogen_structure, pathogen_names, and data
#'   of class `EpiStrainDynamics.pathogen_structure`
#' @family pathogen_structure
#' @export
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#' @srrstats {G2.7, G2.8} data argument accepts standard tabular forms, which
#'   are then post-processed internally to ensure subsequent functions receive
#'   standard inputs
#'
#' @examples
#' subtyped(
#'   data = influenza,
#'   case_timeseries = 'ili',
#'   time = 'week',
#'   influenzaA_unsubtyped_timeseries = 'inf_A',
#'   influenzaA_subtyped_timeseries = c('inf_H3N2', 'inf_H1N1'),
#'   other_pathogen_timeseries = c('inf_B', 'other')
#' )
#'
subtyped <- function (data,
                      case_timeseries,
                      time,
                      influenzaA_unsubtyped_timeseries,
                      influenzaA_subtyped_timeseries,
                      other_pathogen_timeseries) {

  #' @srrstats {G5.8c, G5.8d} edge cases produce errors
  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.1a} validation for vector inputs
  #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling for all numeric
  #'   inputs
  check_column_exists(data, time, "time")
  check_column_exists(data, case_timeseries, "case_timeseries")
  check_column_numeric(data, case_timeseries, "case_timeseries")
  check_column_exists(data, influenzaA_unsubtyped_timeseries,
                      "influenzaA_unsubtyped_timeseries")
  check_column_numeric(data, influenzaA_unsubtyped_timeseries,
                       "influenzaA_unsubtyped_timeseries")

  for (col in influenzaA_subtyped_timeseries) {
    check_column_exists(data, col, "influenzaA_subtyped_timeseries")
    check_column_numeric(data, col, "influenzaA_subtyped_timeseries")
  }

  for (col in other_pathogen_timeseries) {
    check_column_exists(data, col, "other_pathogen_timeseries")
    check_column_numeric(data, col, "other_pathogen_timeseries")
  }

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  all_cols <- c(time, case_timeseries, influenzaA_unsubtyped_timeseries,
                influenzaA_subtyped_timeseries, other_pathogen_timeseries)
  validated_tsbl <- create_validated_tsibble(data, all_cols, time)

  # Create pathogen names
  pathogen_names <- c(influenzaA_subtyped_timeseries, other_pathogen_timeseries)

  # Create matrices directly from validated tsibble
  component_pathogens <- t(as.matrix(
    validated_tsbl[, c(influenzaA_unsubtyped_timeseries, other_pathogen_timeseries)]
  ))

  influenzaA_subtyped <- t(as.matrix(validated_tsbl[, influenzaA_subtyped_timeseries]))

  #' @srrstats {G5.3} data objects are returns with no missing values. this
  #'   is validated through conversion to tsibble in `create_validated_tsibble`
  #' @srrstats {G2.10} extraction of single columns happens after conversion
  #'   to consistent tabular form
  model_inputs <- list(
    pathogen_structure = 'subtyped',
    pathogen_names = pathogen_names,
    data = list(
      case_timeseries = validated_tsbl[[case_timeseries]],
      time = validated_tsbl[[time]],
      component_pathogens = component_pathogens,
      influenzaA_subtyped = influenzaA_subtyped
    )
  )

  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}
