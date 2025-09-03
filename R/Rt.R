#' Generic Method for Rt Analysis
#'
#' S3 generic for computing reproduction numbers from fitted models
#'
#' @param fitted_model Fitted model object with class `EpiStrainDynamics.fit`
#' @param tau_max Integer maximum generation interval in days (default: 7)
#' @param gi_dist Function that returns generation interval probability for given day
#' @param ... Additional arguments passed to metrics calculation
#' @return named list of class `EpiStrainDynamics.metric` containing a dataframe
#'  of the calculated metric outcome (`$measure`), the fit object (`$fit`), and the
#'  constructed model object (`$constructed_model`). The `measure` data frame
#'  contains the median of the epidemiological quantity (`y`), the 50% credible
#'  interval of the quantity (`lb_50` & `ub_50`), the 95% credible interval
#'  (`lb_95` & `ub_95`), the proportion greater than a defined threshold value
#'  (`prop`), the pathogen name (`pathogen`), and the time label (`time`).
#' @export
#'
#' @examples
#' \dontrun{
#'   mod <- construct_model(
#'     method = random_walk(),
#'     pathogen_structure = single(
#'       case_timeseries = sarscov2$cases,
#'       time = sarscov2$date))
#'
#'   fit <- fit_model(mod)
#'
#'   rt <- Rt(mod, tau_max = 7, gi_dist = function(x) 4*x*exp(-2*x))
#' }
#'
Rt <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  validate_class_inherits(fitted_model, 'EpiStrainDynamics.fit')
  UseMethod("Rt")
}

#' @rdname Rt
#' @export
Rt.ps <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_multi_pathogen(fitted_model, tau_max, 'Rt',
                                threshold = 1, use_splines = TRUE,
                                tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', 'EpiStrainDynamics.metric', class(out))
  out
}

#' @rdname Rt
#' @export
Rt.rw <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_multi_pathogen(fitted_model, tau_max, 'Rt',
                                threshold = 1, use_splines = FALSE,
                                tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', 'EpiStrainDynamics.metric', class(out))
  out
}

#' @rdname Rt
#' @export
Rt.ps_single <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_single_pathogen(fitted_model, tau_max, 'Rt',
                                 threshold = 1, use_splines = TRUE,
                                 tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', 'EpiStrainDynamics.metric', class(out))
  out
}

#' @rdname Rt
#' @export
Rt.rw_single <- function(fitted_model, tau_max = 7, gi_dist, ...) {
  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))
  out <- compute_single_pathogen(fitted_model, tau_max, 'Rt',
                                 threshold = 1, use_splines = FALSE,
                                 tau_max = tau_max, gi_dist = gi_dist, g_a = g_a)
  class(out) <- c('Rt', 'EpiStrainDynamics.metric', class(out))
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

  R_denom <- matrix(0, nrow = length(a[, 1]), ncol = 1)
  for(k in 0:(tau_max - 1)) {
    R_denom <- R_denom + exp(a[, time_idx - k]) * gi_dist(k)
  }

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

  R_denom <- matrix(0, nrow = dim(a)[1], ncol = 1)
  for(k in 0:(tau_max - 1)) {
    R_denom <- R_denom + exp(a[, pathogen_idx, time_idx - k]) * gi_dist(k)
  }
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

  R_denom <- matrix(0, nrow = dim(a)[1], ncol = 1)
  for(k in 0:(tau_max - 1)) {
    R_denom <- R_denom + rowSums(exp(a[, , time_idx - k])) * gi_dist(k)
  }
  rowSums(exp(a[, , time_idx])) / (R_denom / g_a)
}
