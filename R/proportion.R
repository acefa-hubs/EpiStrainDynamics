#' Generic Method for proportion Analysis
#'
#' S3 generic for computing proportions from fitted models
#'
#' @param fitted_model Fitted model object with appropriate class
#' @param numerator_combination Named pathogens or subtypes to be included in proportion numerator
#' @param denominator_combination Named pathogens or sutypes to be included in proportion denominator
#' @param ... Additional arguments passed to methods
#' @return Data frame with proportion analysis results
#' @export
proportion <- function(fitted_model,
                       numerator_combination = NULL,
                       denominator_combination = NULL, ...) {
  UseMethod("proportion")
}

#' @rdname proportion
#' @export
proportion.ps <- function(fitted_model,
                          numerator_combination = NULL,
                          denominator_combination = NULL, ...) {

  if (is.null(numerator_combination)) {
    numerator_idx_all <- length(unique(fitted_model$constructed_model$pathogen_names))
    measure <- do.call(rbind, lapply(1:numerator_idx_all,
                                     compute_proportions, fitted_model,
                                     threshold = 0, use_splines = TRUE,
                                     denominator_idx = 1:numerator_idx_all))
  }
  else {
    numerator_idx <- match(numerator_combination, fitted_model$constructed_model$pathogen_names)
    denominator_idx <- match(denominator_combination, fitted_model$constructed_model$pathogen_names)

    measure <- compute_proportions(numerator_idx, fitted_model,
                                   threshold = 0, use_splines = TRUE,
                                   denominator_idx = denominator_idx)
  }

  match_pathogen_name <- function(x) unique(fitted_model$constructed_model$pathogen_names)[x]

  measure <- measure |>
    dplyr::mutate(dplyr::across(dplyr::matches("pathogen"), match_pathogen_name)) |>
    dplyr::rowwise() |>
    dplyr::mutate(pathogen = paste(dplyr::c_across(dplyr::matches('pathogen')), collapse = ', ')) |>
    dplyr::ungroup()

  out <- list(measure = measure,
              fit = fitted_model$fit,
              constructed_model = fitted_model$constructed_model
  )
  class(out) <- c('proportion', class(out))
  out

}

#' @rdname proportion
#' @export
proportion.rw <- function(fitted_model,
                          numerator_combination = NULL,
                          denominator_combination = NULL, ...) {

  if (is.null(numerator_combination)) {
    numerator_idx_all <- length(unique(fitted_model$constructed_model$pathogen_names))
    measure <- do.call(rbind, lapply(1:numerator_idx_all,
                                     compute_proportions, fitted_model,
                                     threshold = 0, use_splines = FALSE,
                                     denominator_idx = 1:numerator_idx_all))
  }
  else {
    numerator_idx <- match(numerator_combination, fitted_model$constructed_model$pathogen_names)
    denominator_idx <- match(denominator_combination, fitted_model$constructed_model$pathogen_names)

    measure <- compute_proportions(numerator_idx, fitted_model,
                                   threshold = 0, use_splines = FALSE,
                                   denominator_idx = denominator_idx)
  }

  match_pathogen_name <- function(x) unique(fitted_model$constructed_model$pathogen_names)[x]

  measure <- measure |>
    dplyr::mutate(dplyr::across(dplyr::matches("pathogen"), match_pathogen_name)) |>
    dplyr::rowwise() |>
    dplyr::mutate(pathogen = paste(dplyr::c_across(dplyr::matches('pathogen')), collapse = ', ')) |>
    dplyr::ungroup()

  out <- list(measure = measure,
              fit = fitted_model$fit,
              constructed_model = fitted_model$constructed_model
  )

  class(out) <- c('proportion', class(out))
  out
}

#' @rdname proportion
#' @export
proportion.ps_single <- function(fitted_model, ...) {
  stop('Proportions only calculated for multi pathogens models.')
}

#' @rdname proportion
#' @export
proportion.rw_single <- function(fitted_model, ...) {
  stop('Proportions only calculated for multi pathogens models.')
}

# =====================
# CALCULATION FUNCTION
# =====================

#' Calculate proportion for Individual Pathogen
#'
#' Computes log-difference in incidence for a specific pathogen
#'
#' @param a Array of log-incidence posterior samples [samples, pathogens, time]
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @param denominator_idx Integer index of the pathogen/s used in denominator of proportion
#' @return Vector of proportion posterior samples
calc_proportion <- function(a, time_idx, pathogen_idx, post, components,
                            denominator_idx) {

  if (length(pathogen_idx) > 1) {
    num <- rowSums(exp(a[, pathogen_idx, time_idx]))
  } else {
    num <- exp(a[, pathogen_idx, time_idx])
  }

  den <- rowSums(exp(a[, denominator_idx, time_idx]))

  num/den
}
