#' Calculate growth rate
#'
#' @param fitted_model fitted model output created with `fit_model`
#'
#' @returns growth rate calculation
#' @export
#'
calculate_growth_rate <- function (fitted_model) {
  gr <- growth_rate(fitted_model)
  return(gr)
}

#' Calculate growth rate
#'
#' @param fitted_model fitted model output created with `fit_model`
#'
#' @returns growth rate calculation
#'
growth_rate <- function (fitted_model) UseMethod("growth_rate")

#' @exportS3Method
growth_rate.rw <- function (fitted_model) {

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
growth_rate.ps <- function (fitted_model) {

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
growth_rate.rw_single <- function (fitted_model) {

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
growth_rate.ps_single <- function (fitted_model) {

  fit <- fitted_model$fit
  days_per_knot <- fitted_model$constructed_model$model_params$days_per_knot
  spline_degree <- fitted_model$constructed_model$model_params$spline_degree
  time <- fitted_model$constructed_model$data$time
  time_labels <- fitted_model$constructed_model$data$time

  B_true <- predict_B_true(time, days_per_knot, spline_degree)

  post <- rstan::extract(fit)

  num_days <- length(time)
  a <- array(data = NA, dim = c(nrow(post$a), num_days))

  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }

  ## multiple ex
  # single_path_gr_internal <- function (obj) {do.call(rbind, lapply(2:num_days, gr_internal, obj, time_labels))}
  #
  # apply(a, 2, single_path_gr_internal)


  df <- do.call(rbind, lapply(2:num_days, gr_internal, obj, time_labels))


  df

}

#' Internal growth rate function
#'
#' @param iter loop id
#' @param obj posterior samples
#' @param time_labels time values
#'
#' @returns row processed for given iter
gr_internal <- function (iter, obj, time_labels) {

  quan <- quantile(obj[, iter] - obj[, (iter - 1)],
                   c(0.5, 0.025, 0.25, 0.75, 0.975))
  prop <- length(obj[, iter][(obj[, iter] - obj[, (iter - 1)]) > 0]) /
    length(obj[, iter] - obj[, (iter - 1)])
  row <- data.frame(time = time_labels[iter],
                    t_step = iter,
                    y = quan[[1]],
                    lb_50 = quan[[3]],
                    ub_50 = quan[[4]],
                    lb_95 = quan[[2]],
                    ub_95 = quan[[5]],
                    prop = prop)
  return(row)
}
