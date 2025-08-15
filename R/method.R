#' Specify random walk method
#'
#' @returns list with method identified as random walk
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
#' @returns list with method and model parameters
#' @export
#'
#' @examples
#' p_spline(spline_degree = 2, days_per_knot = 5)
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
