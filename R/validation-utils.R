#' Validate that input is a positive whole number
#'
#' @param value The value to validate
#' @param arg_name Name of the argument being validated (for error messages)
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
#'
validate_positive_whole_number <- function(value, arg_name) {

  # Check that value is numeric
  if (!is.numeric(value)) {
    cli::cli_abort("Argument {arg_name} must be numeric")
  }

  # Check that value is a single value (length = 1)
  if (length(value) != 1) {
    cli::cli_abort("Argument {arg_name} must be a single value")
  }

  #' @srrstats {G2.16} error to check for undefined values
  # Check that value is finite (not NA, NULL, Inf, -Inf, NaN)
  if (!is.finite(value)) {
    cli::cli_abort("Argument {arg_name} must be a finite number")
  }

  # Check that value is a whole number
  if (value != as.integer(value)) {
    cli::cli_abort("Argument {arg_name} must be a whole number")
  }

  # Check that value is positive
  if (value <= 0) {
    cli::cli_abort("Argument {arg_name} must be a positive number")
  }

  invisible(NULL)
}

#' Validate object class inheritance
#'
#' Validates that an object inherits from one or more specified classes.
#' For multiple classes, the object must inherit from ALL specified classes.
#'
#' @param obj The object to validate
#' @param class_names Character vector of one or more class names to check
#' @param require_all Logical; if TRUE (default), object must inherit from ALL
#'   specified classes. If FALSE, object must inherit from at least ONE class.
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {G5.8, G5.8a} checks for zero length
#'
#' @return Invisible NULL if validation passes, otherwise throws an error
validate_class_inherits <- function(obj, class_names, require_all = TRUE) {

  # Ensure class_names is a character vector
  if (length(class_names) == 0) {
    cli::cli_abort("{class_names} must be a non-empty character vector")
  }

  # Check inheritance for each class
  inheritance_results <- sapply(class_names, function(cls) inherits(obj, cls))

  if (require_all) {
    # Must inherit from ALL classes
    if (!all(inheritance_results)) {
      missing_classes <- class_names[!inheritance_results]
      if (length(class_names) == 1) {
        obj_class <- class(obj)
        cli::cli_abort("Input must be of class {class_names} but got class: {obj_class}")
      } else {
        cli::cli_abort(c("Input must inherit from all classes: {class_names}",
                         "Missing classes: {missing_classes}",
                         "Actual classes: {obj_class}"))
      }
    }
  } else {
    # Must inherit from at least ONE class
    if (!any(inheritance_results)) {
      cli::cli_abort(c("Input must inherit from at least one of: {class_names}",
                       "Actual classes: {obj_class}"))
    }
  }

  invisible(NULL)
}

#' Validate input vectors have matching lengths
#'
#' @param ... Named vectors to check for matching lengths
#' @param .error_call The calling function for better error messages
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
validate_matching_lengths <- function(..., .error_call = rlang::caller_env()) {

  vectors <- list(...)
  vector_names <- names(vectors)
  lengths <- lengths(vectors)

  if (length(unique(lengths)) > 1) {
    length_info <- paste(paste(vector_names, lengths, sep = " (length "), ")", collapse = ", ")
    cli::cli_abort(c("All input vectors must have the same length.",
                     "Found: {length_info}"))
  }
}

#' Validate that a list of vectors all have the same length
#'
#' @param vector_list Named list of vectors
#' @param list_name Name of the list parameter for error messages
#' @param reference_length Expected length (optional)
#' @param .error_call The calling function for better error messages
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {G5.8, G5.8a} checks for zero length
#'
validate_list_vector <- function(vector_list, list_name,
                                 reference_length = NULL,
                                 .error_call = rlang::caller_env()) {

  if (!is.list(vector_list)) {
    cli::cli_abort("{.arg {list_name}} must be a list.")
  }

  if (length(vector_list) == 0) {
    cli::cli_abort("{list_name} cannot be empty")
  }

  lengths <- lengths(vector_list)
  vector_names <- names(vector_list)

  if (is.null(vector_names)) {
    cli::cli_abort("{.arg {list_name}} must have named elements.")
  }

  if (any(vector_names == "" | is.na(vector_names))) {
    cli::cli_abort("All elements in {.arg {list_name}} must have non-empty names.")
  }

  if (length(unique(lengths)) > 1) {
    length_info <- paste(paste(vector_names, lengths, sep = " (length "), ")", collapse = ", ")
    cli::cli_abort(c("All vectors in {list_name} must have the same length.",
                     "Found: {length_info}"))
  }

  if (!is.null(reference_length) && lengths[1] != reference_length) {
    length_obj <- lengths[1]
    cli::cli_abort("Vectors in {list_name} must have length {reference_length} but found length {length_obj}")
  }

  invisible(TRUE)
}

#' Validate input of tau priors
#'
#' @param tau_priors priors object of class `EpiStrainDynamics.priors`
#' @param structure pathogen_structure object, containing information on
#'  the number of pathogens
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {G5.8, G5.8a} checks for zero length
#'
validate_tau_priors <- function(tau_priors, structure) {

  # If tau_priors is NULL, return NULL (valid default)
  if (is.null(tau_priors)) {
    return(NULL)
  }

  # Check if tau_priors has the correct class
  if (!inherits(tau_priors, "EpiStrainDynamics.prior")) {
    cli::cli_abort("tau_priors must be created using the priors() function or be NULL")
  }

  # Check if tau_priors has required components
  if (!all(c("mean", "sd") %in% names(tau_priors))) {
    cli::cli_abort("tau_priors must contain both 'mean' and 'sd' components")
  }

  # Determine expected dimension based on cov_structure
  if (structure$model_params$cov_structure == "0") {
    expected_dim <- 1
  } else if (structure$model_params$cov_structure == "1") {
    if (is.null(structure$pathogen_names)) {
      cli::cli_abort("structure$pathogen_names is required when cov_structure is '1'")
    }
    expected_dim <- length(structure$pathogen_names)
  } else if (structure$model_params$cov_structure == "2") {
    cli::cli_alert("user-supplied priors for tau are not currently supported for correlated covariance structures")
  }

  # Validate and adjust dimensions for mean
  mean_dim <- length(tau_priors$mean)
  if (mean_dim == 1 && expected_dim > 1) {
    # Repeat single value to match expected dimension
    tau_priors$mean <- rep(tau_priors$mean, expected_dim)
  } else if (mean_dim != expected_dim && mean_dim > 1) {
    cli::cli_abort(
      "Length of 'mean' in tau_priors ({mean_dim}) does not match expected dimension ({expected_dim}) for cov_structure '{.var {structure$cov_structure}}'")
  }

  # Validate and adjust dimensions for sd
  sd_dim <- length(tau_priors$sd)
  if (sd_dim == 1 && expected_dim > 1) {
    # Repeat single value to match expected dimension
    tau_priors$sd <- rep(tau_priors$sd, expected_dim)
  } else if (sd_dim != expected_dim && sd_dim > 1) {
    cli::cli_abort(
      "Length of 'sd' in tau_priors ({sd_dim}) does not match expected dimension ({expected_dim}) for cov_structure {.var {structure$cov_structure}}"
    )
  }

  # Return the validated (and potentially adjusted) tau_priors
  return(tau_priors)
}

#' Validate input of phi priors
#'
#' @param phi_priors priors object of class `EpiStrainDynamics.priors`
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {G5.8, G5.8a} checks for zero length
#'
validate_phi_priors <- function(phi_priors) {
  # If phi_priors is NULL, return NULL (valid default)
  if (is.null(phi_priors)) {
    return(NULL)
  }

  # Check if phi_priors has the correct class
  if (!inherits(phi_priors, "EpiStrainDynamics.prior")) {
    cli::cli_abort("phi_priors must be created using the priors() function or be NULL")
  }

  # Check if phi_priors has required components
  if (!all(c("mean", "sd") %in% names(phi_priors))) {
    cli::cli_abort("phi_priors must contain both 'mean' and 'sd' components")
  }

  # Check that both mean and sd are scalar (length 1)
  if (length(phi_priors$mean) != 1) {
    cli::cli_abort("phi_priors$mean must be a single value")
  }

  if (length(phi_priors$sd) != 1) {
    cli::cli_abort("phi_priors$sd must be a single value")
  }

  return(phi_priors)
}
