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
                       adapt_delta = 0.9,
                       multi_cores = TRUE,
                       verbose = TRUE,
                       seed = NULL, ...) {

  # validate inputs
  validate_class_inherits(constructed_model, 'EpiStrainDynamics.model')
  validate_n_chain(n_chain)
  validate_n_iter(n_iter)
  validate_n_warmup(n_warmup)
  validate_thin(thin)
  validate_adapt_delta(adapt_delta)
  validate_seed(seed)
  validate_verbose(verbose)
  validate_multi_cores(multi_cores)

  UseMethod("fit_model")
}

#' @rdname fit_model
#' @export
fit_model.rw_subtyped <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.9,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

  fit_object <- rstan::sampling(
    stanmodels$rw_subtyped,
    data = constructed_model$standata,
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

  # add chain names to init list
  # names(inits) <- paste0("chain", seq(1, n_chain, 1))

  #' @srrstats {BS5.0} function returns initial values used in computation
  #' @srrstats {BS5.5} the `model` return object is of class `stanfit`, which
  #'   includes information about convergence
  out <- list(fit = fit_object,
                      # inits = inits,
                      constructed_model = constructed_model)

  class(out) <- c('rw', 'EpiStrainDynamics.fit', class(out))
  return(out)
}

#' @rdname fit_model
#' @export
fit_model.ps_subtyped <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.9,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

  fit_object <- rstan::sampling(
    stanmodels$ps_subtyped,
    data = constructed_model$standata,
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

#' @rdname fit_model
#' @export
fit_model.rw_multiple <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.9,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

  fit_object <- rstan::sampling(
    stanmodels$rw_multiple,
    data = constructed_model$standata,
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

#' @rdname fit_model
#' @export
fit_model.ps_multiple <- function (constructed_model,
                                   n_chain = 4,
                                   n_iter = 2000,
                                   n_warmup = floor(n_iter/2),
                                   thin = 1,
                                   adapt_delta = 0.9,
                                   multi_cores = TRUE,
                                   verbose = TRUE,
                                   seed = NULL, ...) {

  fit_object <- rstan::sampling(
    stanmodels$ps_multiple,
    data = constructed_model$standata,
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

#' @rdname fit_model
#' @export
fit_model.rw_single <- function (constructed_model,
                                 n_chain = 4,
                                 n_iter = 2000,
                                 n_warmup = floor(n_iter/2),
                                 thin = 1,
                                 adapt_delta = 0.9,
                                 multi_cores = TRUE,
                                 verbose = TRUE,
                                 seed = NULL, ...) {

  fit_object <- rstan::sampling(
    stanmodels$rw_single,
    data = constructed_model$standata,
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

#' @rdname fit_model
#' @export
fit_model.ps_single <- function (constructed_model,
                                 n_chain = 4,
                                 n_iter = 2000,
                                 n_warmup = floor(n_iter/2),
                                 thin = 1,
                                 adapt_delta = 0.9,
                                 multi_cores = TRUE,
                                 verbose = TRUE,
                                 seed = NULL, ...) {

  fit_object <- rstan::sampling(
    stanmodels$ps_single,
    data = constructed_model$standata,
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

