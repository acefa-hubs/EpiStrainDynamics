#' Generic Method for Rt Analysis
#'
#' S3 generic for computing reproduction numbers from fitted models
#'
#' @param fitted_model Fitted model object with appropriate class
#' @param tau_max Integer maximum generation interval in days (default: 7)
#' @param gi_dist Function that returns generation interval probability for given day
#' @param ... Additional arguments passed to methods
#' @return Data frame with Rt analysis results
#' @export
Rt <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  UseMethod("Rt")
}

#' @rdname Rt
#' @export
Rt.ps <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_multi_pathogen(fitted_model, tau_max, 'Rt',
                                threshold = 1, use_splines = TRUE,
                                tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', class(out))
  out
}

#' @rdname Rt
#' @export
Rt.rw <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_multi_pathogen(fitted_model, tau_max, 'Rt',
                                threshold = 1, use_splines = FALSE,
                                tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', class(out))
  out
}

#' @rdname Rt
#' @export
Rt.ps_single <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_single_pathogen(fitted_model, tau_max, 'Rt',
                                 threshold = 1, use_splines = TRUE,
                                 tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', class(out))
  out
}

#' @rdname Rt
#' @export
Rt.rw_single <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_single_pathogen(fitted_model, tau_max, 'Rt',
                                 threshold = 1, use_splines = FALSE,
                                 tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', class(out))
  out
}

# =====================
# CALCULATION FUNCTIONS
# =====================

#' Calculate Rt for Single Pathogen Model
#'
#' Computes reproduction number for single pathogen models
#'
#' @param a Array of log-incidence posterior samples [samples, time]
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @param tau_max Integer maximum generation interval (days)
#' @param gi_dist Function returning generation interval probability for given day
#' @param g_a Numeric normalization constant (sum of generation interval)
#' @return Vector of Rt posterior samples
calc_rt_single <- function(a, time_idx, pathogen_idx, post, components,
                           tau_max, gi_dist, g_a) {
  R_denom <- calc_rt_denominator(a, time_idx, tau_max, gi_dist, FALSE)
  exp(a[, time_idx]) / (R_denom / g_a)
}

#' Calculate Rt for Individual Pathogen
#'
#' Computes reproduction number for a specific pathogen using generation interval convolution
#'
#' @param a Array of log-incidence posterior samples [samples, pathogens, time]
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @param tau_max Integer maximum generation interval (days)
#' @param gi_dist Function returning generation interval probability for given day
#' @param g_a Numeric normalization constant (sum of generation interval)
#' @return Vector of Rt posterior samples
calc_rt_individual <- function(a, time_idx, pathogen_idx, post, components,
                               tau_max, gi_dist, g_a) {
  R_denom <- calc_rt_denominator(a, time_idx, tau_max, gi_dist, TRUE, pathogen_idx)
  exp(a[, pathogen_idx, time_idx]) / (R_denom / g_a)
}

#' Calculate Total Rt Across All Pathogens
#'
#' Computes total reproduction number across all pathogens
#'
#' @param a Array of log-incidence posterior samples [samples, pathogens, time]
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @param tau_max Integer maximum generation interval (days)
#' @param gi_dist Function returning generation interval probability for given day
#' @param g_a Numeric normalization constant (sum of generation interval)
#' @return Vector of total Rt posterior samples
calc_rt_total <- function(a, time_idx, pathogen_idx, post, components,
                          tau_max, gi_dist, g_a) {
  R_denom <- calc_rt_denominator(a, time_idx, tau_max, gi_dist, TRUE, NULL)
  rowSums(exp(a[, , time_idx])) / (R_denom / g_a)
}
