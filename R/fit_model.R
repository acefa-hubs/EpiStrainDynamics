#' Method fit model
#'
#' @param constructed_model prepared model object from construct_model()
#' @param iter number of iterations
#' @param warmup number of warmup iterations
#' @param chains number of chains
#'
#' @returns fit model
#' @export
#'
fit_model <- function (constructed_model,
                       iter = 2000,
                       warmup = 1000,
                       chains = 3) {


  # valid_class_inherits(constructed_model, 'EpiStrainDynamics.')

  UseMethod("fit_model")
}

# fit stan model
#' @exportS3Method
fit_model.rw_subtyped <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names

  standata <- list(num_data = length(cases),
                   num_path = length(pathogen_names),
                   Y = cases,
                   P1 = constructed_model$data$component_pathogens,
                   P2 = constructed_model$data$influenzaA_subtyped,
                   week_effect = constructed_model$model_params$week_effect,
                   DOW = constructed_model$model_params$DOW,
                   cov_structure = constructed_model$model_params$cov_structure,
                   noise_structure = constructed_model$model_params$noise_structure)

  fit_object <- rstan::sampling(stanmodels$rw_subtyped,
                                data = standata,
                                iter = iter,
                                warmup = warmup,
                                chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('rw', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @exportS3Method
fit_model.ps_subtyped <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names
  time_seq <- constructed_model$data$time_seq
  spline_degree <- constructed_model$model_params$spline_degree
  knots <- constructed_model$model_params$knots

  standata <- list(num_data = length(cases),
                   num_knots = length(knots),
                   num_path = length(pathogen_names),
                   knots = knots,
                   spline_degree = spline_degree,
                   Y = cases,
                   P1 = constructed_model$data$component_pathogens,
                   P2 = constructed_model$data$influenzaA_subtyped,
                   X = time_seq,
                   week_effect = constructed_model$model_params$week_effect,
                   DOW = constructed_model$model_params$DOW,
                   cov_structure = constructed_model$model_params$cov_structure,
                   noise_structure = constructed_model$model_params$noise_structure)

  fit_object <- rstan::sampling(stanmodels$ps_subtyped,
                                data = standata,
                                iter = iter,
                                warmup = warmup,
                                chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('ps', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @exportS3Method
fit_model.rw_multiple <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names

  standata <- list(num_data = length(cases),
                   num_path = length(pathogen_names),
                   Y = cases,
                   P = constructed_model$data$component_pathogens,
                   week_effect = constructed_model$model_params$week_effect,
                   DOW = constructed_model$model_params$DOW,
                   cov_structure = constructed_model$model_params$cov_structure,
                   noise_structure = constructed_model$model_params$noise_structure)

  fit_object <- rstan::sampling(stanmodels$rw_multiple,
                                data = standata,
                                iter = iter,
                                warmup = warmup,
                                chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('rw', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @exportS3Method
fit_model.ps_multiple <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names
  time_seq <- constructed_model$data$time_seq
  spline_degree <- constructed_model$model_params$spline_degree
  knots <- constructed_model$model_params$knots

  standata <- list(num_data = length(cases),
                   num_knots = length(knots),
                   num_path = length(pathogen_names),
                   knots = knots,
                   spline_degree = spline_degree,
                   Y = cases,
                   P = constructed_model$data$component_pathogens,
                   X = time_seq,
                   week_effect = constructed_model$model_params$week_effect,
                   DOW = constructed_model$model_params$DOW,
                   cov_structure = constructed_model$model_params$cov_structure,
                   noise_structure = constructed_model$model_params$noise_structure)

  fit_object <- rstan::sampling(stanmodels$ps_multiple,
                                data = standata,
                                iter = iter,
                                warmup = warmup,
                                chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('ps', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @exportS3Method
fit_model.rw_single <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries

  standata <- list(num_data = length(cases),
                   Y = cases,
                   week_effect = constructed_model$model_params$week_effect,
                   DOW = constructed_model$model_params$DOW)

  fit_object <- rstan::sampling(stanmodels$rw_single,
                                data = standata,
                                iter = iter,
                                warmup = warmup,
                                chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('rw_single', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @exportS3Method
fit_model.ps_single <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  time_seq <- constructed_model$data$time_seq
  spline_degree <- constructed_model$model_params$spline_degree
  knots <- constructed_model$model_params$knots

  standata <- list(num_data = length(cases),
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree = spline_degree,
                   Y = cases,
                   X = time_seq,
                   week_effect = constructed_model$model_params$week_effect,
                   DOW = constructed_model$model_params$DOW)

  fit_object <- rstan::sampling(stanmodels$ps_single,
                                data = standata,
                                iter = iter,
                                warmup = warmup,
                                chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('ps_single', 'EpiStrainDynamics.fit', class(out))

  return(out)
}
