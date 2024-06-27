#' Function for getting knot locations
#'
#' @param X days
#' @param days_per_knot days per knot
#' @param spline_degree spline degree
#'
#' @return location of knots
#' @export
#'
get_knots <- function (X, days_per_knot, spline_degree = 3) {

  X <- as.numeric(X)

  num_knots <- ceiling((max(X) - min(X)) / days_per_knot)

  first_knot <- min(X) - spline_degree * days_per_knot
  final_knot <- first_knot + days_per_knot * num_knots +
    2 * spline_degree * days_per_knot

  knots <- seq(first_knot, final_knot, by = days_per_knot)

}
