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

#' Validate prior parameters
#'
#' @param mean The mean parameter value
#' @param sd The sd parameter value
#' @return List with mean and sd, or NULL if both missing
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {G5.8, G5.8a} checks for zero length
#' @srrstats {BS2.3, BS2.4, BS2.5} checks priors for length, appropriate
#'  values, and distributionally appropriate
#'
validate_priors <- function(mean, sd) {

  # Check if both are provided or both are missing
  mean_provided <- !missing(mean) && !is.null(mean)
  sd_provided <- !missing(sd) && !is.null(sd)

  # If one is provided, both must be provided
  if (mean_provided != sd_provided) {
    cli::cli_abort("If specifying priors, both mean and sd must be provided")
  }

  # If neither provided, return NULL
  if (!mean_provided && !sd_provided) {
    return(NULL)
  }

  # Validate parameters when both are provided
  if (!is.numeric(mean) || !is.numeric(sd)) {
    cli::cli_abort("Both {mean} and {sd} must be numeric")
  }

  #' @srrstats {G2.16} error to check for undefined values
  if (any(is.na(mean)) || any(is.na(sd))) {
    cli::cli_abort("{mean} and {sd} cannot contain NA values")
  }

  if (any(mean < 0) || any(sd <= 0)) {
    cli::cli_abort("All {mean} and {sd} values must be positive")
  }

  # Create and return priors object
  priors <- list('mean' = mean, 'sd' = sd)
  class(priors) <- c('EpiStrainDynamics.prior', class(priors))
  return(priors)
}

#' Validate Smoothing Structure Object
#'
#' @param smoothing_obj An object created by \code{smoothing_structure()}
#' @param pathogen_names Character vector of pathogen names (required for "independent" structure)
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {G5.8, G5.8a} checks for zero length
validate_smoothing_structure <- function(smoothing_obj, pathogen_names = NULL) {

  # Check if object has correct class
  if (!inherits(smoothing_obj, "EpiStrainDynamics.smoothing")) {
    stop("smoothing_params must be created using the smoothing_structure() function")
  }

  # For independent structure, validate and adjust dimensions
  if (smoothing_obj$smoothing_type == "independent") {
    if (is.null(pathogen_names)) {
      stop("pathogen_names is required for 'independent' smoothing structure")
    }

    expected_dim <- length(pathogen_names)

    if (!is.null(smoothing_obj$tau_priors)) {
      # Validate and adjust dimensions for mean
      mean_dim <- length(smoothing_obj$tau_priors$mean)
      if (mean_dim == 1 && expected_dim > 1) {
        smoothing_obj$tau_priors$mean <- rep(smoothing_obj$tau_priors$mean, expected_dim)
      } else if (mean_dim != expected_dim && mean_dim > 1) {
        stop(sprintf(
          "Length of tau_mean (%d) does not match number of pathogens (%d)",
          mean_dim, expected_dim
        ))
      }

      # Validate and adjust dimensions for sd
      sd_dim <- length(smoothing_obj$tau_priors$sd)
      if (sd_dim == 1 && expected_dim > 1) {
        smoothing_obj$tau_priors$sd <- rep(smoothing_obj$tau_priors$sd, expected_dim)
      } else if (sd_dim != expected_dim && sd_dim > 1) {
        stop(sprintf(
          "Length of tau_sd (%d) does not match number of pathogens (%d)",
          sd_dim, expected_dim
        ))
      }
    }
  }

  return(smoothing_obj)
}


# Helper: Validate column exists in data
check_column_exists <- function(data, col_name, arg_name) {
  if (!col_name %in% names(data)) {
    stop(sprintf("Column '%s' specified in %s not found in data", col_name, arg_name))
  }
}

# Helper: Validate column is numeric
check_column_numeric <- function(data, col_name, arg_name) {
  if (!is.numeric(data[[col_name]])) {
    stop(sprintf("Column '%s' specified in %s must be numeric, but is %s",
                 col_name, arg_name, class(data[[col_name]])[1]))
  }
}

#' Check for missing data (NA values) in specified columns
#'
#' @param data dataframe to check
#' @param columns character vector of column names to check for missing data
#' @param context character string describing the context (e.g., "case_timeseries")
#'
#' @returns NULL (invisibly) if no missing data, otherwise throws an error
#' @keywords internal
#'
#' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
#' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
#'
check_missing_data <- function(data, columns, context) {
  for (col in columns) {
    if (any(is.na(data[[col]]))) {
      na_count <- sum(is.na(data[[col]]))
      na_indices <- which(is.na(data[[col]]))

      # Show first few indices if there are many
      if (na_count > 5) {
        indices_msg <- paste0(
          "rows ",
          paste(head(na_indices, 5), collapse = ", "),
          ", ... (", na_count, " total)"
        )
      } else {
        indices_msg <- paste0(
          "row(s) ",
          paste(na_indices, collapse = ", ")
        )
      }

      cli::cli_abort(
        "Missing values (NA) found in column {.val {col}} (specified as {context}) at {indices_msg}. All data columns must have complete values (no NAs) before creating time series."
      )
    }
  }
  invisible(NULL)
}

# Helper: Create and validate tsibble from dataframe
create_validated_tsibble <- function(data, columns, time_col) {
  # Subset to only the columns we need
  temp_df <- data[, columns, drop = FALSE]

  # Check for NA values before tsibble conversion
  check_missing_data(temp_df, columns, "input data")

  # Create tsibble
  temp_tsbl <- tryCatch({
    tsibble::tsibble(temp_df, index = !!rlang::sym(time_col))
  }, error = function(e) {
    cli::cli_abort("Error creating tsibble: {e$message}")
  })

  # Check for gaps
  if (tsibble::has_gaps(temp_tsbl)$.gaps) {
    cli::cli_abort("Time series has gaps. All time points must be present with no missing time periods.")
  }

  # Check regularity
  if (!tsibble::is_regular(temp_tsbl)) {
    cli::cli_abort("Time series is irregular. The model requires regularly spaced time intervals.")
  }

  # Check ordering
  if (!tsibble::is_ordered(temp_tsbl)) {
    cli::cli_abort("Time series is not ordered. Please sort the data in chronological order.")
  }

  # Check for empty data
  if (NROW(temp_tsbl) == 0) {
    cli::cli_abort("There is no data to model. Please provide a dataset with at least one observation.")
  }

  return(temp_tsbl)
}
