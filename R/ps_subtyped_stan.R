#' PS Subtyped
#'
#' @export
#'
#' @param num_data number of data points
#' @param num_knots number of knots
#' @param num_path number of pathogens
#' @param knots the sequence of knots
#' @param spline_degree the degree of spline (is equal to order - 1)
#' @param Y outcome variable, eg daily number of cases
#' @param P1 daily number of lab tests positive for influenza A (must be 1st
#'  column) and all other pathogens
#' @param P2 daily number of influenza A subtypes eg H3N2, H1N1
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
ps_subtyped_stan <- function (num_data, num_knots, num_path,
                              knots, spline_degree,
                              Y, P1, P2, X,
                              week_effect, DOW,
                              cov_structure = c(0, 1, 2),
                              noise_structure = c(0, 1),
                              ...) {

  standata <- list(num_data = num_data,
                   num_knots = num_knots,
                   num_path = num_path,
                   knots = knots,
                   spline_degree = spline_degree,
                   Y = Y,
                   P1 = P1,
                   P2 = P2,
                   X = X,
                   week_effect = week_effect,
                   DOW = DOW,
                   cov_structure = cov_structure,
                   noise_structure = noise_structure)

  out <- rstan::sampling(stanmodels$ps_subtyped, data = standata, ...)
  return(out)
}
