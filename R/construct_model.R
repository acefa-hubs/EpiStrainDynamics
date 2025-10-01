#' Construct model
#'
#' @param method either `random_walk()` or `p_spline()`
#' @param pathogen_structure either `single()`, `multiple()`, or `subtyped()`
#' @param smoothing_params argument is optional and defines the structure of the
#'   smoothing terms including optionally setting the smoothing prior tau.
#'   Created with `smoothing_structure()`. NULL option defaults to 'shared'
#'   smoothing structure and default priors.
#' @param dispersion_params argument is optional and defines priors for the
#'   overdispersion parameter of the negative binomial likelihood for
#'   the case timeseries. Created using `dispersion_structure()`. NULL option
#'   uses default priors for phi.
#' @param pathogen_noise logical whether individual pathogen counts have
#'   additional gamma-distributed noise. Default is FALSE. Models with `single`
#'   pathogen structure will be set to FALSE.
#' @param dow_effect logical whether to incorporate a day of week model.
#'
#' @returns a list containing the data, the model parameters, and pathogen
#'  names of class `EpiStrainDynamics.model`
#' @export
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {BS1.1} Code examples of how to enter both simple and more complex
#'   data are here in the `construct_model()` documentation. Descriptive text
#'   can be found in the `README.md`
#' @srrstats {BS1.2, BS1.2c} Code examples of how to specify prior distributions.
#'   Descriptive text can be found in the `README.md`
#' @srrstats {BS2.15} checks that data has been input in correct form, and
#'   provides an informative error message if not
#'
#' @examples
#' # Basic usage
#' method_obj <- random_walk()
#' pathogen_obj <- single(
#'   case_timeseries = sarscov2$cases,
#'   time = sarscov2$date,
#'   pathogen_name = "SARS-CoV-2"
#' )
#' mod <- construct_model(
#'   method = method_obj,
#'   pathogen_structure = pathogen_obj,
#'   dow_effect = TRUE
#' )
#'
#' # Specifying priors on tau and phi are optional
#' mod <- construct_model(
#'
#'   method = p_spline(),
#'
#'   pathogen_structure = multiple(
#'     case_timeseries = sarscov2$cases,
#'     time = sarscov2$date,
#'
#'     component_pathogen_timeseries = list(
#'       alpha = sarscov2$alpha,
#'       delta = sarscov2$delta,
#'       omicron = sarscov2$omicron,
#'       other = sarscov2$other)
#'    ),
#'
#'    smoothing_params = smoothing_structure(
#'       'independent', tau_mean = c(0, 0.1, 0.3, 0), tau_sd = rep(1, times = 4)),
#'    dispersion_params = dispersion_structure(phi_mean = 0, phi_sd = 1),
#'    pathogen_noise = FALSE,
#'    dow_effect = TRUE
#' )
#'
#'
construct_model <- function(method,
                            pathogen_structure,
                            smoothing_params = NULL,
                            dispersion_params = NULL,
                            pathogen_noise = FALSE,
                            dow_effect = FALSE) {

  #' @srrstats {G2.1, G2.2, G5.8, G5.8b} assertions on types of inputs
  validate_class_inherits(
    method, 'EpiStrainDynamics.method'
  )
  validate_class_inherits(
    pathogen_structure, 'EpiStrainDynamics.pathogen_structure'
  )
  smoothing_params <- validate_smoothing_structure(
    smoothing_params, pathogen_structure$pathogen_names)
  if (!is.null(dispersion_params)) {
    validate_class_inherits(dispersion_params, "EpiStrainDynamics.dispersion")
  }

  if (!is.logical(pathogen_noise) || length(pathogen_noise) != 1) {
    cli::cli_abort(
      "You provided {.var {pathogen_noise}} but `pathogen_noise` must be a single
      logical value (TRUE or FALSE)"
    )
  }
  if (!is.logical(dow_effect) || length(dow_effect) != 1) {
    cli::cli_abort(
      "You provided {.var {dow_effect}} but `dow_effect` must be a single
      logical value (TRUE or FALSE)"
    )
  }

  # Extract model type
  model_type <- get_model_type(method$method,
                               pathogen_structure$pathogen_structure)

  # Create time sequence
  time_seq <- seq_len(length(pathogen_structure$data$case_timeseries))

  # Set up day-of-week effects
  week_effect <- if (dow_effect) 7L else 1L
  DOW <- (time_seq - 1L) %% week_effect + 1L

  cases <- pathogen_structure$data$case_timeseries
  pathogen_names <- pathogen_structure$pathogen_names
  component_pathogens <- pathogen_structure$data$component_pathogens %||% NULL
  influenzaA_subtyped <- pathogen_structure$data$influenzaA_subtyped %||% NULL
  cov_structure <- get_cov_structure[smoothing_params$smoothing_type] %||% NULL
  noise_structure <- as.numeric(pathogen_noise) %||% NULL
  spline_degree <- method$model_params$spline_degree %||% NULL

  # priors
  tau_priors_provided <- ifelse(is.null(tau_priors), 1, 2)
  phi_priors_provided <- ifelse(is.null(phi_priors), 1, 2)
  tau_mean <- smoothing_params$tau_priors$mean %||% NULL
  tau_sd <- smoothing_params$tau_priors$sd %||% NULL
  phi_mean <- dispersion_params$phi_priors$mean %||% NULL
  phi_sd <- dispersion_params$phi_priors$sd %||% NULL

  standata <- list(num_data = length(cases),
                   Y = cases,
                   week_effect = week_effect,
                   DOW = DOW,
                   tau_priors_provided = tau_priors_provided,
                   tau_mean = tau_mean,
                   tau_sd = tau_sd,
                   phi_priors_provided = phi_priors_provided,
                   phi_mean = phi_mean,
                   phi_sd = phi_sd)

  if (pathogen_structure$pathogen_structure == 'subtyped') {
    standata <- c(standata,
                  list(num_path = length(pathogen_names),
                       cov_structure = cov_structure,
                       noise_structure = noise_structure,
                       P1 = component_pathogens,
                       P2 = influenzaA_subtyped)
    )
  }

  if (pathogen_structure$pathogen_structure == 'multiple') {
    standata <- c(standata,
                  list(num_path = length(pathogen_names),
                       cov_structure = cov_structure,
                       noise_structure = noise_structure,
                       P = component_pathogens)
    )
  }

  if (method$method == 'p-spline') {
    knots <- get_knots(
      time_seq,
      days_per_knot = method$model_params$days_per_knot,
      spline_degree = method$model_params$spline_degree
    )
    standata <- c(standata,
                  list(num_knots = length(knots),
                       knots = knots,
                       spline_degree = spline_degree,
                       X = time_seq)
    )
  }

  # Prepare data list
  data <- c(
    list(time_seq = time_seq),
    pathogen_structure$data
  )

  # return priors
  # tau_priors = smoothing_params$tau_priors,
  # phi_priors = dispersion_params$phi_priors,

  # Construct final model input list
  model_input <- list(
    data = data,
    standata = standata,
    pathogen_names = pathogen_structure$pathogen_names,
    dow_effect = dow_effect
  )

  class(model_input) <- c(model_type, "EpiStrainDynamics.model",
                          class(model_input))
  return(model_input)
}

#' Extract model type from method and pathogen structure
#'
#' @param method_name Character string: method name ('random-walk' or 'p-spline')
#' @param pathogen_type Character string: pathogen structure type
#'
#' @noRd
#' @returns Character string: model type
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
get_model_type <- function(method_name, pathogen_type) {

  # Input validation
  valid_methods <- c('random-walk', 'p-spline')
  valid_pathogens <- c('single', 'multiple', 'subtyped')

  if (!method_name %in% valid_methods) {
    valid_methods_pasted <- paste(valid_methods, collapse = ", ")
    cli::cli_abort(
      "Unknown `method`: {.var {method_name}}.
      Valid methods include {.var {valid_methods_pasted}}"
    )
  }

  if (!pathogen_type %in% valid_pathogens) {
    valid_pathogens_pasted <- paste(valid_pathogens, collapse = ", ")
    cli::cli_abort(
      "Unknown `pathogen_type`: {.var {pathogen_type}}.
      Valid pathogen types include {.var {valid_pathogens_pasted}}"
    )
  }

  # Lookup model type
  method_abbrev <- switch(method_name,
                          'p-spline' = 'ps',
                          'random-walk' = 'rw')

  model_type <- paste(method_abbrev, pathogen_type, sep = '_')

  return(model_type)
}

#' Function for getting knot locations
#'
#' @param X Numeric vector of time points
#' @param days_per_knot Number of days between knots (must be positive)
#' @param spline_degree Polynomial degree of spline (must be positive)
#'
#' @noRd
#' @return Numeric vector of knot locations
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
get_knots <- function(X, days_per_knot = 3, spline_degree = 3) {

  # Input validation
  if (!is.numeric(X) || length(X) == 0) {
    cli::cli_abort(
      "Argument `X` must be a non-empty numeric vector"
    )
  }

  validate_positive_whole_number(days_per_knot, "days_per_knot")
  validate_positive_whole_number(spline_degree, "spline_degree")

  X <- as.numeric(X)
  min_X <- min(X, na.rm = TRUE)
  max_X <- max(X, na.rm = TRUE)

  if (is.na(min_X) || is.na(max_X)) {
    cli::cli_abort(
      "All values in `X` are missing"
    )
  }

  # Calculate knot parameters
  num_knots <- ceiling((max_X - min_X) / days_per_knot)
  first_knot <- min_X - spline_degree * days_per_knot
  final_knot <- first_knot + days_per_knot * (num_knots + 2L * spline_degree)

  knots <- seq(first_knot, final_knot, by = days_per_knot)

  return(knots)
}

#' Get covariance structure from assigned smoothing structure
#'
#' @param smoothing_structure either `shared` (all pathogens have the same
#'   smoothing structure; \code{tau[1]}), `independent` (each pathogen with
#'   independent smoothing structure; \code{tau[number of pathogens]}), or
#'   `correlated` (smoothing structure is correlated among pathogens
#'   \code{Sigma[number of pathogens, number of pathogens]}). Case-insensitive.
#'
#' @returns numeric value for covariance structure needed by stan models
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
get_cov_structure <- function(smoothing_structure = c('shared',
                                                      'independent',
                                                      'correlated')) {

  # Convert to numeric codes
  cov_structure <- switch(smoothing_structure,
                          'shared' = 0,
                          'independent' = 1,
                          'correlated' = 2,                            {
                            cli::cli_abort("Invalid option provided: '{smoothing_structure}'.
                                             Please choose 'shared', 'independent', or 'correlated'.")
                          }
  )

  return(cov_structure)
}
