#' Specify random walk method
#'
#' @returns list with method identified as random walk of class
#'   `EpiStrainDynamics.method`
#' @export
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
#' @param spline_degree polynomial degree of the individual spline segments
#'   used to construct the overall curve (must be a positive whole number)
#' @param days_per_knot number of days for each knot (must be a positive whole number)
#'
#' @returns list with method and model parameters of class
#'   `EpiStrainDynamics.method`
#' @export
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

  # Validate input parameters
  validate_positive_whole_number(spline_degree, "spline_degree")
  validate_positive_whole_number(days_per_knot, "days_per_knot")

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
