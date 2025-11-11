#' Specify random walk method
#'
#' Random walk is one of two optional Bayesian smoothing prior methods that can
#' be selected and used in the model definition with `EpiStrainDynamics`. The
#' random walk is a stochastic process that describes a path of a series of
#' random steps on a mathematical space, and each next step's direction only
#' depends on the current position, not the previous path. No additional
#' arguments must be supplied to define the random walk method.
#'
#' @returns list with method identified as random walk of class
#'   `EpiStrainDynamics.method`
#' @family method
#' @export
#'
#' @srrstats {G1.3} method defined clearly
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#'
#' @examples random_walk()
random_walk <- function () {
  model_inputs <- list(
    method = 'random-walk'
  )

  class(model_inputs) <- 'EpiStrainDynamics.method'
  return(model_inputs)
}

#' Specify p_spline method
#'
#' Penalised splines is one of two optional Bayesian smoothing prior methods
#' that can be selected and used in the model definition with
#' `EpiStrainDynamics`. The main benefit of selecting the penalised spline
#' method over random walks is that it can capture dynamical effects (fine
#' enough temporal resolution) while not being too computationally expensive.
#'
#' @param spline_degree polynomial degree of the individual spline segments
#'   used to construct the overall curve (must be a positive whole number)
#' @param days_per_knot number of days for each knot (must be a positive whole number)
#'
#' @returns list with method and model parameters of class
#'   `EpiStrainDynamics.method`
#' @family method
#' @export
#'
#' @srrstats {G1.3} method defined clearly
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#'
#' @examples
#' # Valid usage
#' p_spline(spline_degree = 2L, days_per_knot = 5L)
#' p_spline(spline_degree = 3, days_per_knot = 7)
#'
#' # These will produce validation errors (as intended):
#' \donttest{
#' # Non-positive values
#' try(p_spline(spline_degree = 0, days_per_knot = 5))
#' try(p_spline(spline_degree = 3, days_per_knot = -1))
#'
#' # Non-whole numbers
#' try(p_spline(spline_degree = 2.5, days_per_knot = 5))
#' try(p_spline(spline_degree = 3, days_per_knot = 4.2))
#'
#' # Non-numeric values
#' try(p_spline(spline_degree = "invalid", days_per_knot = 5))
#' }
p_spline <- function (spline_degree = 3,
                      days_per_knot = 3) {

  #' @srrstats {G2.1} assertions on types of inputs
  validate_positive_whole_number(spline_degree, "spline_degree")
  validate_positive_whole_number(days_per_knot, "days_per_knot")

  #' @srrstats {G2.4a} explicit conversion to `integer` via `as.integer()`
  model_inputs <- list(
    method = 'p-spline',
    model_params = list(
      spline_degree = as.integer(spline_degree),
      days_per_knot = as.integer(days_per_knot)
    )
  )
  class(model_inputs) <- 'EpiStrainDynamics.method'
  return(model_inputs)
}
