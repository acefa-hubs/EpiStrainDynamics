#' Create Prior Specifications for Model Parameters
#'
#' This function creates a standardized prior object for use with
#' `EpiStrainDynamics` model parameters. The prior object contains mean and
#' standard deviation values that define the prior distribution for model
#' parameters.
#'
#' @param mean Numeric vector specifying the prior mean(s). Must be provided.
#' @param sd Numeric vector specifying the prior standard deviation(s). Must be
#'   provided and all values must be positive.
#'
#' @return An object of class \code{EpiStrainDynamics.prior} containing:
#'   \item{mean}{The specified prior mean(s)}
#'   \item{sd}{The specified prior standard deviation(s)}
#'
#' @details
#' Both \code{mean} and \code{sd} must be specified when calling this function.
#' The \code{sd} values must be positive (> 0). No missing values (NA) are allowed
#' in either parameter.
#'
#' The dimensions of \code{mean} and \code{sd} should match the requirements of
#' the specific model parameter being specified. For example:
#' \itemize{
#'   \item For \code{phi_priors}: both \code{mean} and \code{sd} must be scalar (length 1)
#'   \item For \code{tau_priors}: dimensions depend on the covariance structure.
#'   When covariance is 'shared' then a single mean and sd can be provided. When
#'   covariance is 'independent' then the user can provide a different prior for
#'   each pathogen or provide a single one that will be used for each.
#' }
#'
#' @examples
#' # Create scalar priors (e.g., for phi parameter)
#' phi_prior <- priors(mean = 2.0, sd = 0.5)
#'
#' # Create vector priors (e.g., for tau parameter with multiple pathogens)
#' tau_prior <- priors(mean = c(0, 0, 0), sd = c(1, 1, 1))
#'
#' @seealso \code{\link{validate_tau_priors}}, \code{\link{validate_phi_priors}}
#'
#' @export
priors <- function(mean, sd) {
  # Check that both mean and sd are provided
  if (missing(mean) || missing(sd)) {
    stop("Both 'mean' and 'sd' must be specified when using priors()")
  }

  # Check that mean and sd are numeric
  if (!is.numeric(mean) || !is.numeric(sd)) {
    stop("Both 'mean' and 'sd' must be numeric")
  }

  # Check for NA values
  if (any(is.na(mean)) || any(is.na(sd))) {
    stop("'mean' and 'sd' cannot contain NA values")
  }

  # Check that sd values are positive
  if (any(sd <= 0)) {
    stop("All 'sd' values must be positive")
  }

  p_list <- list('mean' = mean, 'sd' = sd)
  class(p_list) <- c('EpiStrainDynamics.prior', class(p_list))
  p_list
}
