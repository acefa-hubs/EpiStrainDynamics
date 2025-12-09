#' Generic Method for proportion Analysis
#'
#' Computes epidemiological proportion, defined as the relative fraction of
#' cases attributable to specific pathogen(s) or strain(s) at a given time point,
#' calculated as the ratio of incidence from selected pathogen(s) to a
#' reference group:
#' \deqn{P_t = \frac{I_{\text{numerator}}(t)}{I_{\text{denominator}}(t)}}
#'
#' Where incidences are derived from the exponential of log-incidence estimates:
#' \deqn{P_t = \frac{\sum_{i \in \text{numerator}} \exp(\log\text{-incidence}_{i,t})}{\sum_{j \in \text{denominator}} \exp(\log\text{-incidence}_{j,t})}}
#'
#' This metric quantifies the **relative contribution** of specific pathogen(s)
#' or strain(s) to the total disease burden. Key characteristics:
#' \itemize{
#'   \item Values between 0 and 1 (\eqn{0 \leq P_t \leq 1}) when denominator
#'    includes numerator components
#'   \item Values can exceed 1 (\eqn{P_t > 1}) when denominator excludes
#'   numerator components
#'   \item Represents the fractional share of cases at each time point
#'   \item Time-varying to capture changing pathogen/strain dynamics
#' }
#'
#' **Flexible combinations**:
#' \itemize{
#'   \item Individual proportions: Each pathogen relative to all pathogens
#'   (default)
#'   \item Custom numerator: Specific pathogen(s) of interest (e.g., variant
#'   of concern)
#'   \item Custom denominator: Either 'all' pathogens or a specified subset
#'   \item Subset comparisons: Compare specific groups (e.g., Alpha vs.
#'   Delta + Omicron)
#' }
#'
#' This metric function can be run directly on the fitted model output.
#'
#' @param fitted_model Fitted model object with class `EpiStrainDynamics.fit`
#'  with `multiple` or `subtyped` pathogen structure.
#' @param numerator_combination Named pathogens or subtypes to be included in
#'  proportion numerator.
#' @param denominator_combination Named pathogens or subtypes to be included in
#'  proportion denominator, or 'all'.
#' @param ... Additional arguments passed to metrics calculation
#' @return named list of class `EpiStrainDynamics.metric` containing a dataframe
#'  of the calculated metric outcome (`$measure`), the fit object (`$fit`), and the
#'  constructed model object (`$constructed_model`). The `measure` data frame
#'  contains the median of the epidemiological quantity (`y`), the 50% credible
#'  interval of the quantity (`lb_50` & `ub_50`), the 95% credible interval
#'  (`lb_95` & `ub_95`), the proportion greater than a defined threshold value
#'  (`prop`), the pathogen name (`pathogen`), and the time label (`time`).
#' @family metrics
#' @export
#'
#' @srrstats {G1.3} metric defined clearly
#' @srrstats {G1.4} uses `Roxygen2` documentation
#'
#' @examples
#' \dontrun{
#'   mod <- construct_model(
#'     method = p_spline(),
#'     pathogen_structure = multiple(
#'       case_timeseries = sarscov2$cases,
#'       time = sarscov2$date,
#'       component_pathogen_timeseries = list(
#'         alpha = sarscov2$alpha,
#'         delta = sarscov2$delta,
#'         omicron = sarscov2$omicron,
#'         other = sarscov2$other))
#'   )
#'
#'   fit <- fit_model(mod)
#'   prop <- proportion(fit)
#'
#'   # or a unique combination, compared to all pathogens
#'   prop2 <- proportion(fit,
#'     numerator_combination = c('alpha', 'delta', 'omicron'),
#'     denominator_combination = 'all'
#'   )
#'
#'   # or a user-specified combination in both numerator and denominator
#'   prop3 <- proportion(fit,
#'     numerator_combination = 'alpha',
#'     denominator_combination = c('alpha', 'delta', 'omicron')
#'   )
#' }
proportion <- function(fitted_model,
                       numerator_combination = NULL,
                       denominator_combination = NULL, ...) {

  validate_class_inherits(fitted_model, 'EpiStrainDynamics.fit')
  validate_class_inherits(fitted_model, c('ps', 'rw'), require_all = FALSE)

  # Validate combination arguments
  pathogen_names <- unique(fitted_model$constructed_model$pathogen_names)
  validate_pathogen_combination(numerator_combination, pathogen_names,
                                "numerator_combination")
  validate_pathogen_combination(denominator_combination, pathogen_names,
                                "denominator_combination")

  UseMethod("proportion")
}

#' @rdname proportion
#' @export
proportion.ps <- function(fitted_model,
                          numerator_combination = NULL,
                          denominator_combination = NULL, ...) {

  if (is.null(numerator_combination)) {
    numerator_idx_all <- length(unique(fitted_model$constructed_model$pathogen_names))

    measure <- do.call(
      rbind, lapply(1:numerator_idx_all,
                    function(x) {
                      compute_single_pathogen(
                        fitted_model, 1, 'proportion',
                        use_splines = TRUE,
                        numerator_idx = x,
                        denominator_idx = 1:numerator_idx_all
                      )
                    }
      ))
  }

  else {
    numerator_idx <- match(numerator_combination, fitted_model$constructed_model$pathogen_names)

    if (identical(denominator_combination, "all")) {
      denominator_idx <- 1:length(unique(fitted_model$constructed_model$pathogen_names))
    }
    else {
      denominator_idx <- match(denominator_combination, fitted_model$constructed_model$pathogen_names)
    }

    measure <- compute_single_pathogen(fitted_model, 1, 'proportion',
                                       use_splines = TRUE,
                                       numerator_idx = numerator_idx,
                                       denominator_idx = denominator_idx)
  }

  match_pathogen_name <- function(x) unique(fitted_model$constructed_model$pathogen_names)[x]

  measure <- measure |>
    dplyr::rowwise() |>
    dplyr::mutate(pathogen = paste(match_pathogen_name(.data$pathogen),
                                   collapse = ', ')) |>
    dplyr::ungroup()

  out <- list(measure = measure,
              fit = fitted_model$fit,
              constructed_model = fitted_model$constructed_model
  )
  class(out) <- c('proportion', 'EpiStrainDynamics.metric', class(out))
  out

}

#' @rdname proportion
#' @export
proportion.rw <- function(fitted_model,
                          numerator_combination = NULL,
                          denominator_combination = NULL, ...) {

  if (is.null(numerator_combination)) {
    numerator_idx_all <- length(unique(fitted_model$constructed_model$pathogen_names))

    measure <- do.call(
      rbind, lapply(1:numerator_idx_all,
                    function(x) {
                      compute_single_pathogen(
                        fitted_model, 1, 'proportion',
                        use_splines = FALSE,
                        numerator_idx = x,
                        denominator_idx = 1:numerator_idx_all
                      )
                    }
      ))
  }

  else {
    numerator_idx <- match(numerator_combination, fitted_model$constructed_model$pathogen_names)

    if (denominator_combination == 'all') {
      denominator_idx <- 1:length(unique(fitted_model$constructed_model$pathogen_names))
    }
    else {
      denominator_idx <- match(denominator_combination, fitted_model$constructed_model$pathogen_names)
    }

    measure <- compute_single_pathogen(fitted_model, 1, 'proportion',
                                       use_splines = FALSE,
                                       numerator_idx = numerator_idx,
                                       denominator_idx = denominator_idx)
  }

  match_pathogen_name <- function(x) {
    unique(fitted_model$constructed_model$pathogen_names)[x]
  }
  measure <- measure |>
    dplyr::mutate(dplyr::across(dplyr::matches("pathogen"), match_pathogen_name)) |>
    dplyr::rowwise() |>
    dplyr::mutate(pathogen = paste(dplyr::c_across(dplyr::matches('pathogen')),
                                   collapse = ', ')) |>
    dplyr::ungroup()

  out <- list(measure = measure,
              fit = fitted_model$fit,
              constructed_model = fitted_model$constructed_model
  )
  class(out) <- c('proportion', 'EpiStrainDynamics.metric', class(out))
  out
}

# =====================
# CALCULATION FUNCTION
# =====================

#' Calculate proportion for Individual Pathogen
#'
#' Computes log-difference in incidence for a specific pathogen
#'
#' @param a Array of log-incidence posterior samples \code{[samples, pathogens, time]}
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @param numerator_idx Integer index of the pathogen/s used in numerator of proportion
#' @param denominator_idx Integer index of the pathogen/s used in denominator of proportion
#'
#' @return Vector of proportion posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_proportion <- function(a, time_idx, pathogen_idx, post, components,
                            numerator_idx, denominator_idx) {

  if (length(pathogen_idx) > 1) {
    num <- rowSums(exp(a[, pathogen_idx, time_idx]))
  } else {
    num <- exp(a[, pathogen_idx, time_idx])
  }

  den <- rowSums(exp(a[, denominator_idx, time_idx]))

  num/den
}
