#' Create Smoothing Structure Specification with Priors
#'
#' This function creates a standardized smoothing structure object that specifies
#' both the smoothing structure and associated priors for EpiStrainDynamics models.
#'
#' @param smoothing_type Character string specifying the smoothing type:
#'   \itemize{
#'     \item \code{"shared"}: All pathogens have the same smoothing parameter
#'       (equivalent to \code{tau[1]}). By default a model with a single pathogen
#'       will have `shared` smoothing type.
#'     \item \code{"independent"}: Independent smoothing per pathogen
#'       (equivalent to \code{tau[number of pathogens]})
#'     \item \code{"correlated"}: Correlated smoothing type (equivalent
#'       to \code{Sigma[number of pathogens, number of pathogens]})
#'   }
#' @param tau_mean Optional numeric vector specifying the prior mean(s) for tau
#'   parameter. Can be provided for `shared` (single value) and `independent`
#'   smoothing types (can provide a single value which will be repeated for
#'   each pathogen or can provide a unique prior for each pathogen). Prior for
#'   tau for `correlated` smoothing type is not currently supported.
#' @param tau_sd Numeric vector specifying the prior standard deviation(s)
#'   for tau parameter. Can be provided for `shared` (single value) and `independent`
#'   smoothing types (can provide a single value which will be repeated for
#'   each pathogen or can provide a unique prior for each pathogen). Prior for
#'   tau for `correlated` smoothing type is not currently supported.
#'
#' @return An object of class `EpiStrainDynamics.smoothing` containing:
#'   \item{smoothing_type}{The specified smoothing structure type}
#'   \item{tau_priors}{Prior specifications for tau}
#'
#' @srrstats {G1.3} clear definitions of smoothing structure types and prior
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS1.2, BS1.2c} Code examples of how to specify prior distributions.
#'   Descriptive text can be found in the `README.md`
#' @srrstats {BS2.2} Prior distributions are pre-processed and validated
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#'
#' @examples
#' # Shared smoothing with scalar priors
#' shared_smooth <- smoothing_structure("shared", tau_mean = 0, tau_sd = 1)
#'
#' # Independent smoothing with vector priors
#' indep_smooth <- smoothing_structure("independent",
#'                                     tau_mean = c(0, 0, 0),
#'                                     tau_sd = c(1, 1, 1))
#'
#' @export
smoothing_structure <- function(smoothing_type,
                                tau_mean = NULL,
                                tau_sd = NULL) {

  #' @srrstats {G2.3, G2.3a, G2.3b} univariate variables must match and be
  #'   case insensitive
  valid_structures <- c("shared", "independent", "correlated")
  smoothing_type <- tolower(rlang::arg_match(
    arg = smoothing_type,
    values = valid_structures
  ))

  # Handle priors based on structure type
  if (smoothing_type == "correlated") {
    # No priors needed for correlated structure
    if (!is.null(tau_mean) || !is.null(tau_sd)) {
      cli::cli_alert("{tau_mean} and {tau_sd} are ignored for 'correlated' smoothing type")
    }
    tau_priors <- NULL
  } else {
    # For shared and independent structures, validate priors
    tau_priors <- validate_priors(tau_mean, tau_sd)
  }

  #' @srrstats {BS5.2} return priors
  # Create the smoothing structure object
  smooth_obj <- list(
    smoothing_type = smoothing_type,
    tau_priors = tau_priors
  )

  class(smooth_obj) <- c('EpiStrainDynamics.smoothing', class(smooth_obj))
  return(smooth_obj)
}

#' Create Dispersion Structure Specification
#'
#' This function creates a dispersion structure object that specifies priors
#' for the overdispersion parameter of the negative binomial likelihood for
#' the case timeseries.
#'
#' @param phi_mean Numeric scalar specifying the prior mean for the negative
#'   binomial dispersion parameter. Must be positive. Optional - if not
#'   provided, no priors will be set.
#' @param phi_sd Numeric scalar specifying the prior standard deviation for the
#'   dispersion parameter. Must be positive. Optional - if not provided, no
#'   priors will be set.
#'
#' @return An object of class `EpiStrainDynamics.dispersion` containing prior
#'   specifications for phi parameter
#'
#' @srrstats {G1.3} clear definitions of dispersion parameters
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS1.2, BS1.2c} Code examples of how to specify prior distributions.
#'   Descriptive text can be found in the `README.md`
#' @srrstats {BS2.2} Prior distributions are pre-processed and validated
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#'
#' @examples
#' # Create dispersion structure
#' disp_struct <- dispersion_structure(phi_mean = 2.0, phi_sd = 0.5)
#'
#' @export
dispersion_structure <- function(phi_mean = NULL, phi_sd = NULL) {

  # Validate dispersion parameters using helper function
  disp_obj <- validate_priors(phi_mean, phi_sd)

  # Additional check for dispersion priors (must be scalar)
  if (length(phi_mean) != 1 || length(phi_sd) != 1) {
    cli::cli_abort("Only one value for {phi_mean} and {phi_sd} expected.")
  }

  #' @srrstats {BS5.2} return priors
  class(disp_obj) <- c('EpiStrainDynamics.dispersion', class(disp_obj))
  return(disp_obj)
}
