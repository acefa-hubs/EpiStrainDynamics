#' Single pathogen structure
#'
#' @param data dataframe containing columns with all relevant data, or a time
#'   series object (ts, xts, zoo, tsibble, etc.)
#' @param case_timeseries Column name containing case counts. Must be numeric
#'   or a \code{units} object from the \pkg{units} package.
#' @param time name of column with time data. Required for non-time-series
#'   input data. Flexible format - can be date, index, or others, accepted
#'   as `index` identifiers in the `tsibble` time format. Optional when
#'   `data` is a time series class object (ts, mts, xts, zoo, zooreg, tsibble)
#'   as the time index will be automatically detected.
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
#' @srrstats {TS1.7} accommodate units data input
#'
#' @examples
#' # Using a data frame
#' single(
#'   data = sarscov2,
#'   case_timeseries = 'cases',
#'   time = 'date'
#' )
#'
#' # Using a time series object (time argument is optional)
#' \dontrun{
#' sarscov2_xts <- xts::xts(sarscov2[, c("cases", "alpha")], order.by = sarscov2$date)
#' single(
#'   data = sarscov2_xts,
#'   case_timeseries = 'cases'
#' )
#' }
#'
single <- function(data,
                   case_timeseries,
                   time = NULL) {

  #' @srrstats {G5.8c, G5.8d} edge cases produce errors
  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.1a} validation for vector inputs
  if(!is.null(time)) check_column_exists(data, time)
  check_column_exists(data, case_timeseries)

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  validated_tsbl <- create_validated_timeseries(data, case_timeseries, time)
  colnames(validated_tsbl)[colnames(validated_tsbl) == case_timeseries] <-
    "case_timeseries"

  #' @srrstats {G5.3} data objects are returns with no missing values. this
  #'   is validated through conversion to tsibble in `create_validated_timeseries`
  #' @srrstats {G2.10} extraction of single columns happens after conversion
  #'   to consistent tabular form
  model_inputs <- list(
    pathogen_structure = 'single',
    pathogen_names = case_timeseries,
    validated_tsbl = validated_tsbl,
    data = list(
      case_timeseries = validated_tsbl$case_timeseries
    )
  )

  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}

#' Multiple pathogen structure
#'
#' @param data dataframe containing columns with all relevant data, or a time
#'   series object (ts, xts, zoo, tsibble, etc.)
#' @param case_timeseries Column name containing case counts. Must be numeric
#'   or a \code{units} object from the \pkg{units} package.
#' @param time name of column with time data. Required for non-time-series
#'   input data. Flexible format - can be date, index, or others, accepted
#'   as `index` identifiers in the `tsibble` time format. Optional when
#'   `data` is a time series class object (ts, mts, xts, zoo, zooreg, tsibble)
#'   as the time index will be automatically detected.
#' @param component_pathogen_timeseries vector of column names with additional
#'   pathogen case count timeseries to model. Must be numeric or a \code{units}
#'   object from the \pkg{units} package.
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
#' @srrstats {TS1.7} accommodate units data input
#'
#' @examples
#' # Using a data frame
#' multiple(
#'   data = sarscov2,
#'   case_timeseries = 'cases',
#'   time = 'date',
#'   component_pathogen_timeseries = c('alpha', 'delta', 'omicron', 'other')
#' )
#'
#' # Using a time series object (time argument is optional)
#' \dontrun{
#' sarscov2_xts <- xts::xts(
#'   sarscov2[, c("cases", "alpha", "delta", "omicron", "other")],
#'   order.by = sarscov2$date
#' )
#' multiple(
#'   data = sarscov2_xts,
#'   case_timeseries = 'cases',
#'   component_pathogen_timeseries = c('alpha', 'delta', 'omicron', 'other')
#' )
#' }
#'
multiple <- function(data,
                     case_timeseries,
                     time = NULL,
                     component_pathogen_timeseries) {

  #' @srrstats {G5.8c, G5.8d} edge cases produce errors
  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.1a} validation for vector inputs
  if(!is.null(time)) check_column_exists(data, time)
  ts_cols <- c(case_timeseries, component_pathogen_timeseries)
  for (col in ts_cols) {
    check_column_exists(data, col)
  }

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  validated_tsbl <- create_validated_timeseries(data, ts_cols, time)
  colnames(validated_tsbl)[colnames(validated_tsbl) == case_timeseries] <-
    "case_timeseries"

  # Create matrix from validated tsibble
  component_pathogens <- t(as.matrix(
    validated_tsbl[, component_pathogen_timeseries])
  )

  #' @srrstats {G5.3} data objects are returned with no missing values.
  #'   this is validated through conversion to tsibble in
  #'   `create_validated_timeseries`
  #' @srrstats {G2.10} extraction of single columns happens after conversion
  #'   to consistent tabular form
  model_inputs <- list(
    pathogen_structure = 'multiple',
    pathogen_names = component_pathogen_timeseries,
    validated_tsbl = validated_tsbl,
    data = list(
      case_timeseries = validated_tsbl$case_timeseries,
      component_pathogens = component_pathogens
    )
  )
  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}

#' Subtyped pathogen structure
#'
#' @param data dataframe containing columns with all relevant data, or a time
#'   series object (ts, xts, zoo, tsibble, etc.)
#' @param case_timeseries Column name containing case counts. Must be numeric
#'   or a \code{units} object from the \pkg{units} package.
#' @param time name of column with time data. Required for non-time-series
#'   input data. Flexible format - can be date, index, or others, accepted
#'   as `index` identifiers in the `tsibble` time format. Optional when
#'   `data` is a time series class object (ts, mts, xts, zoo, zooreg, tsibble)
#'   as the time index will be automatically detected.
#' @param influenzaA_unsubtyped_timeseries vector of column names with additional
#'   unsubtyped influenzaA case count timeseries. Must be numeric or a
#'   \code{units} object from the \pkg{units} package.
#' @param influenzaA_subtyped_timeseries vector of column names with additional
#'   subtyped influenzaA case count timeseries. Must be numeric or a
#'   \code{units} object from the \pkg{units} package.
#' @param other_pathogen_timeseries vector of column names with additional
#'   pathogen case count timeseries to model. Must be numeric or a \code{units}
#'   object from the \pkg{units} package.
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
#' @srrstats {TS1.7} accommodate units data input
#'
#' @examples
#' # Using a data frame
#' subtyped(
#'   data = influenza,
#'   case_timeseries = 'ili',
#'   time = 'week',
#'   influenzaA_unsubtyped_timeseries = 'inf_A',
#'   influenzaA_subtyped_timeseries = c('inf_H3N2', 'inf_H1N1'),
#'   other_pathogen_timeseries = c('inf_B', 'other')
#' )
#'
#' # Using a time series object (time argument is optional)
#' \dontrun{
#' influenza_xts <- xts::xts(
#'   influenza[, c("ili", "inf_A", "inf_H3N2", "inf_H1N1", "inf_B", "other")],
#'   order.by = influenza$week
#' )
#' subtyped(
#'   data = influenza_xts,
#'   case_timeseries = 'ili',
#'   influenzaA_unsubtyped_timeseries = 'inf_A',
#'   influenzaA_subtyped_timeseries = c('inf_H3N2', 'inf_H1N1'),
#'   other_pathogen_timeseries = c('inf_B', 'other')
#' )
#' }
#'
subtyped <- function(data,
                     case_timeseries,
                     time = NULL,
                     influenzaA_unsubtyped_timeseries,
                     influenzaA_subtyped_timeseries,
                     other_pathogen_timeseries) {

  #' @srrstats {G5.8c, G5.8d} edge cases produce errors
  #' @srrstats {G2.1, G5.8, G5.8b} assertions on types of inputs
  #' @srrstats {G2.1a} validation for vector inputs
  if(!is.null(time)) check_column_exists(data, time)
  ts_cols <- c(case_timeseries, influenzaA_unsubtyped_timeseries,
                influenzaA_subtyped_timeseries, other_pathogen_timeseries)
  for (col in ts_cols) {
    check_column_exists(data, col)
  }

  #' @srrstats {G2.0, G2.0a} validate matching lengths for input data
  #' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
  #' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
  validated_tsbl <- create_validated_timeseries(data, ts_cols, time)
  colnames(validated_tsbl)[colnames(validated_tsbl) == case_timeseries] <-
    "case_timeseries"

  # Create pathogen names
  pathogen_names <- c(influenzaA_subtyped_timeseries, other_pathogen_timeseries)

  # Create matrices directly from validated tsibble
  component_pathogens <- t(as.matrix(
    validated_tsbl[, c(influenzaA_unsubtyped_timeseries,
                       other_pathogen_timeseries)]
  ))

  influenzaA_subtyped <- t(as.matrix(
    validated_tsbl[, influenzaA_subtyped_timeseries])
  )

  #' @srrstats {G5.3} data objects are returns with no missing values. this
  #'   is validated through conversion to tsibble in `create_validated_timeseries`
  #' @srrstats {G2.10} extraction of single columns happens after conversion
  #'   to consistent tabular form
  model_inputs <- list(
    pathogen_structure = 'subtyped',
    pathogen_names = pathogen_names,
    validated_tsbl = validated_tsbl,
    data = list(
      case_timeseries = validated_tsbl$case_timeseries,
      component_pathogens = component_pathogens,
      influenzaA_subtyped = influenzaA_subtyped
    )
  )

  class(model_inputs) <- c('EpiStrainDynamics.pathogen_structure',
                           class(model_inputs))
  model_inputs
}
