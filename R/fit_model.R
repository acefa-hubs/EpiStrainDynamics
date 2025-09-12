#' Generic Method for fitting model
#'
#' S3 generic for fitted models from constructed model object
#' describe how to set up inits.
#'
#' @srrstats {BS2.7, BS2.11} Starting values can be specified with the
#'   `initial_values` argument.
#'
#' @param constructed_model prepared model object of class
#'  `EpiStrainDynamics.model`
#' @param n_chain number of MCMC chains, defaults to 4
#' @param n_iter A positive integer specifying the number of iterations for each chain, default value is  2000
#' @param n_warmup A positive integer specifying the number of warmup iterations, default value is half the number of iterations
#' @param thin A positive integer specifying the period for saving samples, default value is 1.
#' @param adapt_delta Numeric value between 0 and 1 indicating target average acceptance probability used in `rstan::sampling`. Default value is 0.8.
#' @param multi_cores A logical value indicating whether to parallelize chains with multiple cores, default is TRUE and uses all available cores - 1.
#' @param verbose Logical value controlling the verbosity of output (i.e., warnings, messages, progress bar), default is TRUE.
#' @param seed A positive integer seed used for random number generation in MCMC. Default is NULL, which means the seed is generated from 1 to the maximum integer supported by R.
#' @param ... additional arguments to `rstan::sampling()`
#'
#' @returns fit model of class `EpiStrainDynamics.fit`
#' @export
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.3, BS1.3} Arguments `n_chain`, `n_warmup`, `n_iter`, `thin`, and
#'   `adapt_delta` control the computational process and are described clearly
#'   in the documentation.
#' @srrstats {BS2.12} Argument `verbose`, which defaults to `TRUE` controls the
#'   verbosity of the stan sampling output.
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#' @srrstats {BS5.0} *Return values should include starting value(s) or seed(s), including values for each sequence where multiple sequences are included*
#' @srrstats {BS5.1} Returned list is of class `EpiStrainDynamics.fit`, and
#'   includes the fit object of class `stanfit` and the pre-specified
#'   constructed model object of class `EpiStrainDynamics.model`. The
#'   model object itself includes the input data.
#'
#'
#' @examples
#' \dontrun{
#'   mod <- construct_model(
#'     method = random_walk(),
#'     pathogen_structure = single(
#'       case_timeseries = sarscov2$cases,
#'       time = sarscov2$date))
#'
#'   fit <- fit_model(mod)
#'
#'   # or specify additional mcmc parameters
#'   fit <- fit_model(
#'     mod, iter = 3000, warmup = 2000, chains = 4
#'   )
#' }
#'
fit_model <- function (constructed_model,
                       n_chain = 4,
                       n_iter = 2000,
                       n_warmup = floor(n_iter/2),
                       thin = 1,
                       adapt_delta = 0.8,
                       multi_cores = TRUE,
                       verbose = TRUE,
                       seed = NULL, ...) {

  validate_class_inherits(constructed_model, 'EpiStrainDynamics.model')
  UseMethod("fit_model")
}

#' @rdname growth_rate
#' @export
fit_model.rw_subtyped <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.8,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

  cases <- constructed_model$data$case_timeseries
  pathogen_names <- constructed_model$pathogen_names

  standata <- list(
    num_data = length(cases),
    num_path = length(pathogen_names),
    Y = cases,
    P1 = constructed_model$data$component_pathogens,
    P2 = constructed_model$data$influenzaA_subtyped,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    cov_structure = constructed_model$model_params$cov_structure,
    noise_structure = constructed_model$model_params$noise_structure
  )

  fit_object <- rstan::sampling(
    stanmodels$rw_subtyped,
    data = standata,
    #' @srrstats {G2.4, G2.4a} explicit conversion to integers
    chains = as.integer(n_chain),
    iter = as.integer(n_iter),
    warmup = as.integer(n_warmup),
    thin = as.integer(thin),
    seed = ifelse(!is.null(seed), as.integer(seed),
                  sample.int(.Machine$integer.max, 1)),
    cores = ifelse(multi_cores == TRUE, parallel::detectCores()-1, 1),
    control = list(adapt_delta = adapt_delta),
    refresh = ifelse(verbose == TRUE, 500, 0),
    ...
  )

  # check for divergent transitions
  div_trans <- sum(lapply(rstan::get_sampler_params(fit_object,
                                                    inc_warmup = FALSE),
                          div_check)[[1]])
  ## print either troubleshooting or visualization tips
  # if (div_trans > 0 && verbose) {
  #   url <- "https://ednajoint.netlify.app/tips#troubleshooting-tips"
  #   message <- "Refer to the eDNAjoint guide for troubleshooting tips: "
  # } else {
  #   url <- "https://ednajoint.netlify.app/tips#visualization-tips"
  #   message <- "Refer to the eDNAjoint guide for visualization tips: "
  # }
  # cat(message, url, "\n")

  # add chain names to init list
  # names(inits) <- paste0("chain", seq(1, n_chain, 1))

  #' @srrstats {BS5.0} function returns initial values used in computation
  #' @srrstats {BS5.5} the `model` return object is of class `stanfit`, which
  #'   includes information about convergence
  result_list <- list(fit = fit_object,
                      # inits = inits,
                      constructed_model = constructed_model)

  class(result_list) <- c('rw', 'EpiStrainDynamics.fit', class(out))

  return(result_list)
}

#' @rdname growth_rate
#' @export
fit_model.ps_subtyped <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.8,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

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

  fit_object <- rstan::sampling(
    stanmodels$ps_subtyped,
    data = standata,
    chains = as.integer(n_chain),
    iter = as.integer(n_iter),
    warmup = as.integer(n_warmup),
    thin = as.integer(thin),
    seed = ifelse(!is.null(seed), as.integer(seed),
                  sample.int(.Machine$integer.max, 1)),
    cores = ifelse(multi_cores == TRUE, parallel::detectCores()-1, 1),
    control = list(adapt_delta = adapt_delta),
    refresh = ifelse(verbose == TRUE, 500, 0),
    ...
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('ps', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @rdname growth_rate
#' @export
fit_model.rw_multiple <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.8,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

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

  fit_object <- rstan::sampling(
    stanmodels$rw_multiple,
    data = standata,
    chains = as.integer(n_chain),
    iter = as.integer(n_iter),
    warmup = as.integer(n_warmup),
    thin = as.integer(thin),
    seed = ifelse(!is.null(seed), as.integer(seed),
                  sample.int(.Machine$integer.max, 1)),
    cores = ifelse(multi_cores == TRUE, parallel::detectCores()-1, 1),
    control = list(adapt_delta = adapt_delta),
    refresh = ifelse(verbose == TRUE, 500, 0),
    ...
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('rw', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @rdname growth_rate
#' @export
fit_model.ps_multiple <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.8,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

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

  fit_object <- rstan::sampling(
    stanmodels$ps_multiple,
    data = standata,
    chains = as.integer(n_chain),
    iter = as.integer(n_iter),
    warmup = as.integer(n_warmup),
    thin = as.integer(thin),
    seed = ifelse(!is.null(seed), as.integer(seed),
                  sample.int(.Machine$integer.max, 1)),
    cores = ifelse(multi_cores == TRUE, parallel::detectCores()-1, 1),
    control = list(adapt_delta = adapt_delta),
    refresh = ifelse(verbose == TRUE, 500, 0),
    ...
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('ps', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @rdname growth_rate
#' @export
fit_model.rw_single <- function (constructed_model,
                                 n_chain = 4,
                                 n_iter = 2000,
                                 n_warmup = floor(n_iter/2),
                                 thin = 1,
                                 adapt_delta = 0.8,
                                 multi_cores = TRUE,
                                 verbose = TRUE,
                                 seed = NULL, ...) {

  cases <- constructed_model$data$case_timeseries

  standata <- list(num_data = length(cases),
                   Y = cases,
                   week_effect = constructed_model$model_params$week_effect,
                   DOW = constructed_model$model_params$DOW)

  fit_object <- rstan::sampling(
    stanmodels$rw_single,
    data = standata,
    chains = as.integer(n_chain),
    iter = as.integer(n_iter),
    warmup = as.integer(n_warmup),
    thin = as.integer(thin),
    seed = ifelse(!is.null(seed), as.integer(seed),
                  sample.int(.Machine$integer.max, 1)),
    cores = ifelse(multi_cores == TRUE, parallel::detectCores()-1, 1),
    control = list(adapt_delta = adapt_delta),
    refresh = ifelse(verbose == TRUE, 500, 0),
    ...
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('rw_single', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' @rdname growth_rate
#' @export
fit_model.ps_single <- function (constructed_model,
                                 n_chain = 4,
                                 n_iter = 2000,
                                 n_warmup = floor(n_iter/2),
                                 thin = 1,
                                 adapt_delta = 0.8,
                                 multi_cores = TRUE,
                                 verbose = TRUE,
                                 seed = NULL, ...) {

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

  fit_object <- rstan::sampling(
    stanmodels$ps_single,
    data = standata,
    chains = as.integer(n_chain),
    iter = as.integer(n_iter),
    warmup = as.integer(n_warmup),
    thin = as.integer(thin),
    seed = ifelse(!is.null(seed), as.integer(seed),
                  sample.int(.Machine$integer.max, 1)),
    cores = ifelse(multi_cores == TRUE, parallel::detectCores()-1, 1),
    control = list(adapt_delta = adapt_delta),
    refresh = ifelse(verbose == TRUE, 500, 0),
    ...
  )

  out <- list(fit = fit_object,
              constructed_model = constructed_model)

  class(out) <- c('ps_single', 'EpiStrainDynamics.fit', class(out))

  return(out)
}

#' Extract Model Components
#'
#' Extracts components from a fitted model object for analysis functions
#'
#' @param model_obj Model object
#'
#' @return List containing extracted components: fit, pathogen_names, num_path,
#'   time_seq, time, num_days, days_per_knot, spline_degree, DOW, week_effect,
#'   dow_effect
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
get_model_components <- function(model_obj) {
  list(
    pathogen_names = fitted_model$constructed_model$pathogen_names %||% NULL,
    num_path = length(fitted_model$constructed_model$pathogen_names %||% 1),
    time_seq = fitted_model$constructed_model$data$time_seq,
    time = fitted_model$constructed_model$data$time,
    num_days = length(fitted_model$constructed_model$data$time),
    days_per_knot = fitted_model$constructed_model$model_params$days_per_knot %||% NULL,
    spline_degree = fitted_model$constructed_model$model_params$spline_degree %||% NULL,
    DOW = fitted_model$constructed_model$model_params$DOW %||% NULL,
    week_effect = fitted_model$constructed_model$model_params$week_effect %||% NULL,
    dow_effect = fitted_model$constructed_model$dow_effect
  )
}
