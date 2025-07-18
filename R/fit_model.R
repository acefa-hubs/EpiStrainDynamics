#' Fit model
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

  # if(!'EpiStrainDynamics.options' %in% class(options))
  #   stop("`options` must be created with model_options()")

  out <- fit(constructed_model, iter, warmup, chains)

  return(out)
}

#' Method fit
#'
#' @param constructed_model prepared model object from construct_model()
#' @param iter number of iterations
#' @param warmup number of warmup iterations
#' @param chains number of chains
#'
#' @returns fit model
#'
fit <- function (constructed_model, iter, warmup, chains) UseMethod("fit")

# fit stan model
#' @exportS3Method
fit.rw_subtyped <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names

  fit_object <- rw_subtyped_stan(
    num_data = length(cases),
    num_path = length(pathogen_names),
    Y = cases,
    P1 = constructed_model$data$component_pathogens,
    P2 = constructed_model$data$influenzaA_subtyped,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    cov_structure = constructed_model$model_params$cov_structure,
    noise_structure = constructed_model$model_params$noise_structure,
    iter = iter,
    warmup = warmup,
    chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
fit.ps_subtyped <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names
  time <- constructed_model$data$time
  spline_degree <- constructed_model$model_params$spline_degree
  knots <- get_knots(
    time,
    days_per_knot = constructed_model$model_params$days_per_knot,
    spline_degree = spline_degree
  )

  fit_object <- ps_subtyped_stan(
    num_data = length(cases),
    num_path = length(pathogen_names),
    num_knots = length(knots),
    knots = knots,
    spline_degree = spline_degree,
    Y = cases,
    X = time,
    P1 = constructed_model$data$component_pathogens,
    P2 = constructed_model$data$influenzaA_subtyped,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    cov_structure = constructed_model$model_params$cov_structure,
    noise_structure = constructed_model$model_params$noise_structure,
    iter = iter,
    warmup = warmup,
    chains = chains)

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
fit.rw_multiple <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names

  fit_object <- rw_multiple_stan(
    num_data = length(cases),
    num_path = length(pathogen_names),
    Y = cases,
    P = constructed_model$data$component_pathogens,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    cov_structure = constructed_model$model_params$cov_structure,
    noise_structure = constructed_model$model_params$noise_structure,
    iter = iter,
    warmup = warmup,
    chains = chains
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
fit.ps_multiple <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names
  time <- constructed_model$data$time
  spline_degree <- constructed_model$model_params$spline_degree
  knots <- get_knots(
    time,
    days_per_knot = constructed_model$model_params$days_per_knot,
    spline_degree = spline_degree
  )

  fit_object <- ps_multiple_stan(
    num_data = length(cases),
    num_path = length(pathogen_names),
    num_knots = length(knots),
    knots = knots,
    spline_degree = spline_degree,
    Y = cases,
    X = time,
    P = constructed_model$data$component_pathogens,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    cov_structure = constructed_model$model_params$cov_structure,
    noise_structure = constructed_model$model_params$noise_structure,
    iter = iter,
    warmup = warmup,
    chains = chains
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
fit.rw_single <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries

  fit_object <- rw_single_stan(
    num_data = length(cases),
    Y = cases,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    iter = iter,
    warmup = warmup,
    chains = chains
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}

#' @exportS3Method
fit.ps_single <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$case_timeseries
  time <- constructed_model$data$time
  spline_degree <- constructed_model$model_params$spline_degree
  knots <- get_knots(
    time,
    days_per_knot = constructed_model$model_params$days_per_knot,
    spline_degree = spline_degree
  )

  fit_object <- ps_single_stan(
    num_data = length(cases),
    num_knots = length(knots),
    knots = knots,
    spline_degree = spline_degree,
    Y = cases,
    X = time,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    iter = iter,
    warmup = warmup,
    chains = chains
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}
