#' Generic Method for Incidence Analysis
#'
#' S3 generic for computing incidence from fitted models
#'
#' @param fitted_model Fitted model object with appropriate class
#' @param ... Additional arguments passed to methods
#' @return Data frame with incidence analysis results
#' @export
incidence <- function(fitted_model, ...) {
  UseMethod("incidence")
}

#' Incidence Analysis for Penalized Spline Multi-Pathogen Models
#'
#' @param fitted_model Fitted model object of class 'ps'
#' @param ... Additional arguments (unused)
#' @return Data frame with incidence results for individual pathogens and total
#' @method incidence ps
#' @export
incidence.ps <- function(fitted_model, ...) {
  compute_multi_pathogen(fitted_model, 1, 'incidence',
                         threshold = 0, use_splines = TRUE)
}

#' Incidence Analysis for Random Walk Multi-Pathogen Models
#'
#' @param fitted_model Fitted model object of class 'rw'
#' @param ... Additional arguments (unused)
#' @return Data frame with incidence results for individual pathogens and total
#' @method incidence rw
#' @export
incidence.rw <- function(fitted_model, ...) {
  compute_multi_pathogen(fitted_model, 1, 'incidence',
                         threshold = 0, use_splines = FALSE)
}

#' Incidence Analysis for Penalized Spline Single Pathogen Models
#'
#' @param fitted_model Fitted model object of class 'ps_single'
#' @param ... Additional arguments (unused)
#' @return Data frame with incidence results
#' @method incidence ps_single
#' @export
incidence.ps_single <- function(fitted_model, ...) {
  compute_single_pathogen(fitted_model, 1, 'incidence',
                          threshold = 0, use_splines = TRUE)
}

#' Incidence Analysis for Random Walk Single Pathogen Models
#'
#' @param fitted_model Fitted model object of class 'rw_single'
#' @param ... Additional arguments (unused)
#' @return Data frame with incidence results
#' @method incidence rw_single
#' @export
incidence.rw_single <- function(fitted_model, ...) {
  compute_single_pathogen(fitted_model, 1, 'incidence',
                          threshold = 0, use_splines = FALSE)
}

# =====================
# CALCULATION FUNCTIONS
# =====================

#' Calculate Incidence for Single Pathogen Model (with DOW adjustment)
#'
#' Computes incidence for single pathogen models, applying day-of-week effects if present
#'
#' @param a Array of log-incidence posterior samples [samples, time]
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object containing day_of_week_simplex if DOW effects used
#' @param components Model components containing DOW information
#' @return Vector of incidence posterior samples
calc_incidence_single <- function(a, time_idx, pathogen_idx, post, components) {
  incidence <- exp(a[, time_idx])

  # Apply DOW effect if required
  if (components$dow_effect) {
    incidence <- incidence * components$week_effect * post$day_of_week_simplex[, components$DOW[time_idx]]
  }

  incidence
}

#' Calculate Incidence for Individual Pathogen (with DOW adjustment)
#'
#' Computes incidence for a specific pathogen, applying day-of-week effects if present
#'
#' @param a Array of log-incidence posterior samples [samples, pathogens, time]
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object containing day_of_week_simplex if DOW effects used
#' @param components Model components containing DOW information
#' @return Vector of incidence posterior samples
calc_incidence_individual <- function(a, time_idx, pathogen_idx, post, components) {
  incidence <- exp(a[, pathogen_idx, time_idx])

  # Apply DOW effect if required
  if (components$dow_effect) {
    incidence <- incidence * components$week_effect *
      post$day_of_week_simplex[, components$DOW[time_idx]]
  }

  incidence
}

#' Calculate Total Incidence Across All Pathogens (with DOW adjustment)
#'
#' Computes total incidence across all pathogens, applying day-of-week effects if present
#'
#' @param a Array of log-incidence posterior samples [samples, pathogens, time]
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object containing day_of_week_simplex if DOW effects used
#' @param components Model components containing DOW information
#' @return Vector of total incidence posterior samples
calc_incidence_total <- function(a, time_idx, pathogen_idx, post, components) {
  incidence <- rowSums(exp(a[, , time_idx]))

  # Apply DOW effect if required
  if (components$dow_effect) {
    incidence <- incidence * components$week_effect *
      post$day_of_week_simplex[, components$DOW[time_idx]]
  }

  incidence
}
