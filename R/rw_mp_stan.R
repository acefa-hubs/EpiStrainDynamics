#' RW MP
#'
#' @export
#'
#' @param num_data number of data points
#' @param num_path number of pathogens
#' @param Y daily number of 'cases'
#' @param P daily number of lab tests positive for each pathogen
#' @param week_effect Number of days in day of week effect? 1=none, 2=weekends?, 7=all days
#' @param DOW integer of day of the week
#' @param cov_structure 0 is tau[1], 1 is tau[num_path], 2 is Sigma[num_path, num_path]
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
rw_mp_stan <- function (num_data, num_path,
                        Y, P, week_effect, DOW,
                        cov_structure = c(0, 1, 2), ...) {

  standata <- list(num_data = num_data,
                   num_path = num_path,
                   Y = Y,
                   P = P,
                   week_effect = week_effect,
                   DOW = DOW,
                   cov_structure = cov_structure)

  out <- rstan::sampling(stanmodels$rw_mp, data = standata, ...)
  return(out)
}
