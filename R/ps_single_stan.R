#' PS single
#'
#' @export
#'
#' @param num_data number of data points
#' @param num_knots num of knots
#' @param knots the sequence of knots
#' @param spline_degree the degree of spline (is equal to order - 1)
#' @param Y
#' @param X
#' @param week_effect Number of days in day of week effect? 1=none, 2=weekends?, 7=all days
#' @param DOW integer of day of the week
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
ps_single_stan <- function (num_data, num_knots, knots,
                            spline_degree, Y, X,
                            week_effect, DOW, ...) {

  standata <- list(num_data = num_data,
                   num_knots = num_knots,
                   knots = knots,
                   spline_degree = spline_degree,
                   Y = Y,
                   X = X,
                   week_effect = week_effect,
                   DOW = DOW)

  out <- rstan::sampling(stanmodels$ps_single, data = standata, ...)
  return(out)
}
