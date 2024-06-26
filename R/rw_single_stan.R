#' RW single
#'
#' @export
#'
#' @param num_data number of data points
#' @param Y daily number of 'cases'
#' @param week_effect number of days in day of week effect? 1=none, 2=weekends?, 7=all days
#' @param DOW integer of day of the week
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
rw_single_stan <- function (num_data, Y, week_effect, DOW, ...) {

  standata <- list(num_data = num_data,
                   Y = Y,
                   week_effect = week_effect,
                   DOW = DOW)

  out <- rstan::sampling(stanmodels$rw_single, data = standata, ...)
  return(out)
}
