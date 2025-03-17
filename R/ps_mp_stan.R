#' PS MP
#'
#' @export
#'
#' @param num_data number of data points
#' @param num_knots number of knots
#' @param num_path number of pathogens
#' @param knots the sequence of knots
#' @param spline_degree the degree of spline (is equal to order - 1)
#' @param Y outcome variable, eg daily number of cases
#' @param P daily number of lab tests positive for each pathogen
#' @param X time data
#' @param week_effect Number of distinct days in day of week effect. 1 = single
#'  effect shared by all days (essentially no DOW), 7 = each day a unique effect
#' @param DOW integer of day of the week
#' @param cov_structure 0 is tau of 1, 1 is tau of
#'  num_path, 2 is Sigma with two dimensions: num_path, num_path.
#' @param noise_structure 0 only includes observation noise (same between
#'  pathogens), 1 includes noise in individual pathogens as well
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
ps_mp_stan <- function (num_data, num_knots, num_path,
                        knots, spline_degree,
                        Y, P, X, week_effect, DOW,
                        cov_structure = c(0, 1, 2),
                        noise_structure = c(0, 1),
                        ...) {

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
                   cov_structure = cov_structure,
                   noise_structure = noise_structure)

  out <- rstan::sampling(stanmodels$ps_mp, data = standata, ...)
  return(out)
}
