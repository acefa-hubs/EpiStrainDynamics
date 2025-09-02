#' Generic Method for Incidence Analysis
#'
#' S3 generic for computing incidence from fitted models
#'
#' @param fitted_model Fitted model object with appropriate class
#' @param dow Logical whether or not to include day-of-week in incidence calc
#' @param ... Additional arguments passed to methods
#' @return Data frame with incidence analysis results
#' @export
incidence <- function(fitted_model, dow, ...) {
  validate_class_inherits(fitted_model, 'EpiStrainDynamics.fit')
  if (dow & !fitted_model$constructed_model$dow_effect) {
    stop("dow effects can't be incorporated into incidence as it was not specified in the model")
  }
  UseMethod("incidence")
}

#' @rdname incidence
#' @export
incidence.ps <- function(fitted_model, dow, ...) {
  out <- compute_multi_pathogen(fitted_model, 1, 'incidence',
                                threshold = 0, use_splines = TRUE, dow)
  class(out) <- c('incidence', class(out))
  out
}

#' @rdname incidence
#' @export
incidence.rw <- function(fitted_model, dow, ...) {
  out <- compute_multi_pathogen(fitted_model, 1, 'incidence',
                                threshold = 0, use_splines = FALSE, dow)
  class(out) <- c('incidence', class(out))
  out
}

#' @rdname incidence
#' @export
incidence.ps_single <- function(fitted_model, dow, ...) {
  out <- compute_single_pathogen(fitted_model, 1, 'incidence',
                                 threshold = 0, use_splines = TRUE, dow)
  class(out) <- c('incidence', class(out))
  out
}

#' @rdname incidence
#' @export
incidence.rw_single <- function(fitted_model, dow, ...) {
  out <- compute_single_pathogen(fitted_model, 1, 'incidence',
                                 threshold = 0, use_splines = FALSE, dow)
  class(out) <- c('incidence', class(out))
  out
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
#' @param dow_effect Logical whether day-of-week effect to be calculated
#' @return Vector of incidence posterior samples
calc_incidence_single <- function(a, time_idx, pathogen_idx,
                                  post, components, dow_effect) {
  incidence <- exp(a[, time_idx])

  # Apply DOW effect if required
  if (dow_effect) {
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
#' @param dow_effect Logical whether day-of-week effect to be calculated
#' @return Vector of incidence posterior samples
calc_incidence_individual <- function(a, time_idx, pathogen_idx,
                                      post, components, dow_effect) {
  incidence <- exp(a[, pathogen_idx, time_idx])

  # Apply DOW effect if required
  if (dow_effect) {
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
#' @param dow_effect Logical whether day-of-week effect to be calculated
#' @return Vector of total incidence posterior samples
calc_incidence_total <- function(a, time_idx, pathogen_idx,
                                 post, components, dow_effect) {
  incidence <- rowSums(exp(a[, , time_idx]))

  # Apply DOW effect if required
  if (dow_effect) {
    incidence <- incidence * components$week_effect *
      post$day_of_week_simplex[, components$DOW[time_idx]]
  }

  incidence
}
