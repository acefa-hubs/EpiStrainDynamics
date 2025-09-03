#' Generic Method for Growth Rate Analysis
#'
#' S3 generic for computing growth rates from fitted models
#'
#' @param fitted_model Fitted model object with class `EpiStrainDynamics.fit`
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
#'   gr <- growth_rate(mod)
#' }
growth_rate <- function(fitted_model, ...) {
  validate_class_inherits(fitted_model, 'EpiStrainDynamics.fit')
  UseMethod("growth_rate")
}

#' @rdname growth_rate
#' @export
growth_rate.ps <- function(fitted_model, ...) {
  out <- compute_multi_pathogen(fitted_model, 2, 'growth_rate',
                                threshold = 0, use_splines = TRUE)
  class(out) <- c('growth_rate', 'EpiStrainDynamics.metric', class(out))
  out
}

#' @rdname growth_rate
#' @export
growth_rate.rw <- function(fitted_model, ...) {
  out <- compute_multi_pathogen(fitted_model, 2, 'growth_rate',
                                threshold = 0, use_splines = FALSE)
  class(out) <- c('growth_rate', 'EpiStrainDynamics.metric', class(out))
  out
}

#' @rdname growth_rate
#' @export
growth_rate.ps_single <- function(fitted_model, ...) {
  out <- compute_single_pathogen(fitted_model, 2, 'growth_rate',
                                 threshold = 0, use_splines = TRUE)
  class(out) <- c('growth_rate', 'EpiStrainDynamics.metric', class(out))
  out
}

#' @rdname growth_rate
#' @export
growth_rate.rw_single <- function(fitted_model, ...) {
  out <- compute_single_pathogen(fitted_model, 2, 'growth_rate',
                                 threshold = 0, use_splines = FALSE)
  class(out) <- c('growth_rate', 'EpiStrainDynamics.metric', class(out))
  out
}

# =====================
# CALCULATION FUNCTIONS
# =====================

#' Calculate Growth Rate for Single Pathogen Model
#'
#' Computes log-difference in incidence for single pathogen models
#'
#' @param a Array of log-incidence posterior samples [samples, time]
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @return Vector of growth rate posterior samples
calc_growth_single <- function(a, time_idx, pathogen_idx, post, components) {
  a[, time_idx] - a[, time_idx - 1]
}

#' Calculate Growth Rate for Individual Pathogen
#'
#' Computes log-difference in incidence for a specific pathogen
#'
#' @param a Array of log-incidence posterior samples [samples, pathogens, time]
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @return Vector of growth rate posterior samples
calc_growth_individual <- function(a, time_idx, pathogen_idx, post, components) {
  a[, pathogen_idx, time_idx] - a[, pathogen_idx, time_idx - 1]
}

#' Calculate Total Growth Rate Across All Pathogens
#'
#' Computes log-difference in total incidence across all pathogens
#'
#' @param a Array of log-incidence posterior samples [samples, pathogens, time]
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @return Vector of total growth rate posterior samples
calc_growth_total <- function(a, time_idx, pathogen_idx, post, components) {
  total_i <- rowSums(exp(a[, , time_idx]))
  total_im1 <- rowSums(exp(a[, , time_idx - 1]))
  log(total_i) - log(total_im1)
}
