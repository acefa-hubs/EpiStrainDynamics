#' Validate Smoothing Structure Object
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {G5.8, G5.8a} checks for zero length
validate_smoothing_structure <- function(smoothing_obj,
                                         pathogen_str,
                                         pathogen_names = NULL) {

  if (!inherits(smoothing_obj, "EpiStrainDynamics.smoothing")) {
    stop("smoothing_params must be created using the smoothing_structure() function")
  }

  smoothing_obj <- validate_smoothing_for_single(smoothing_obj, pathogen_str)
  smoothing_obj <- validate_smoothing_independent(smoothing_obj, pathogen_names)
  smoothing_obj <- validate_smoothing_shared(smoothing_obj)
  smoothing_obj <- validate_smoothing_correlated(smoothing_obj)

  return(smoothing_obj)
}

#' @noRd
validate_smoothing_for_single <- function(smoothing_obj, pathogen_str) {
  if (pathogen_str != "single") return(smoothing_obj)

  if (smoothing_obj$smoothing_type != "shared") {
    smoothing_obj$smoothing_type <- "shared"
    cli::cli_alert('smoothing_type can only be "shared" for "single" pathogen_type')
  }
  if (is.null(smoothing_obj$tau_priors$mean)) {
    smoothing_obj$tau_priors$mean <- 0.0
    smoothing_obj$tau_priors$sd <- 1.0
  }
  return(smoothing_obj)
}

#' @noRd
validate_smoothing_independent <- function(smoothing_obj, pathogen_names) {
  if (smoothing_obj$smoothing_type != "independent") return(smoothing_obj)

  if (is.null(pathogen_names)) {
    stop("pathogen_names is required for 'independent' smoothing structure")
  }
  expected_dim <- length(pathogen_names)

  if (smoothing_obj$priors_provided == 1) {
    smoothing_obj$tau_priors$mean <- rep(0.0, expected_dim)
    smoothing_obj$tau_priors$sd   <- rep(1.0, expected_dim)
  } else {
    smoothing_obj$tau_priors$mean <- validate_smoothing_dim(
      smoothing_obj$tau_priors$mean, expected_dim, "tau_mean"
    )
    smoothing_obj$tau_priors$sd <- validate_smoothing_dim(
      smoothing_obj$tau_priors$sd, expected_dim, "tau_sd"
    )
  }
  return(smoothing_obj)
}

#' @noRd
validate_smoothing_dim <- function(values, expected_dim, name) {
  actual_dim <- length(values)
  if (actual_dim == 1 && expected_dim > 1) {
    return(rep(values, expected_dim))
  }
  if (actual_dim != expected_dim && actual_dim > 1) {
    stop(sprintf(
      "Length of %s (%d) does not match number of pathogens (%d)",
      name, actual_dim, expected_dim
    ))
  }
  return(values)
}

#' @noRd
validate_smoothing_shared <- function(smoothing_obj) {
  if (smoothing_obj$smoothing_type != "shared") return(smoothing_obj)
  if (smoothing_obj$priors_provided != 1) return(smoothing_obj)

  if (is.null(smoothing_obj$tau_priors$mean) ||
      length(smoothing_obj$tau_priors$mean) == 0) {
    smoothing_obj$tau_priors$mean <- array(0.0, dim = 1)
    smoothing_obj$tau_priors$sd   <- array(1.0, dim = 1)
  }
  return(smoothing_obj)
}

#' @noRd
validate_smoothing_correlated <- function(smoothing_obj) {
  if (smoothing_obj$smoothing_type != "correlated") return(smoothing_obj)
  if (smoothing_obj$priors_provided != 1) return(smoothing_obj)

  if (is.null(smoothing_obj$tau_priors$mean) ||
      length(smoothing_obj$tau_priors$mean) == 0) {
    smoothing_obj$tau_priors$mean <- numeric(0)
    smoothing_obj$tau_priors$sd   <- numeric(0)
  }
  return(smoothing_obj)
}
