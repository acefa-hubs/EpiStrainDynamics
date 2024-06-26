#' PS MP
#'
#' @export
#'
#' @param num_data number of data points
#' @param num_knots num of knots
#' @param num_path
#' @param knots the sequence of knots
#' @param spline_degree the degree of spline (is equal to order - 1)
#' @param Y
#' @param P
#' @param X
#' @param week_effect Number of days in day of week effect? 1=none, 2=weekends?, 7=all days
#' @param DOW[num_data] integer of day of the week
#' @param cov_structure 0 is tau[1], 1 is tau[num_path], 2 is Sigma[num_path, num_path]
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
ps_mp_stan <- function (num_data, num_knots, num_path,
                        knots, spline_degree,
                        Y, P, X, week_effect, DOW,
                        cov_structure = c(0, 1, 2), ...) {

  standata <- list(num_data = num_data,
                   num_knots = num_knots,
                   num_path = num_path,
                   knots = knots,
                   spline_degree = spline_degree,
                   Y = Y,
                   P = P,
                   X = X,
                   week_effect = week_effect,
                   DOW = DOW,
                   cov_structure = cov_structure)

  out <- rstan::sampling(stanmodels$ps_mp, data = standata, ...)
  return(out)
}
