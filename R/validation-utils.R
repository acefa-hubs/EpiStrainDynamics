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
  inheritance_results <- vapply(class_names, function(cls) inherits(obj, cls),
                                FUN.VALUE = logical(1))

  obj_class <- class(obj)

  if (require_all) {
    # Must inherit from ALL classes
    if (!all(inheritance_results)) {
      missing_classes <- class_names[!inheritance_results]
      if (length(class_names) == 1) {
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

  # If both are provided, validate them
  if (mean_provided && sd_provided) {
    # Validate parameters
    if (!is.numeric(mean) || !is.numeric(sd)) {
      cli::cli_abort("Both {.var mean} and {.var sd} must be numeric")
    }

    #' @srrstats {G2.16} error to check for undefined values
    if (any(is.na(mean)) || any(is.na(sd))) {
      cli::cli_abort("{.var mean} and {.var sd} cannot contain NA values")
    }

    # Check that mean and sd have the same length
    if (length(mean) != length(sd)) {
      cli::cli_abort(
        "{.var mean} and {.var sd} must have the same length",
        "i" = "{.var mean} has length {length(mean)}, {.var sd} has length {length(sd)}"
      )
    }

    if (any(mean < 0) || any(sd <= 0)) {
      cli::cli_abort("All {.var mean} and {.var sd} values must be positive")
    }
  }

  # Always create priors list structure (even if both are NULL)
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
validate_smoothing_structure <- function(smoothing_obj, pathogen_names = NULL,
                                         pathogen_str) {
  # Check if object has correct class
  if (!inherits(smoothing_obj, "EpiStrainDynamics.smoothing")) {
    stop("smoothing_params must be created using the smoothing_structure() function")
  }

  if (pathogen_str == 'single' & smoothing_obj$smoothing_type != 'shared') {
    smoothing_obj$smoothing_type <- 'shared'
    cli::cli_alert('smoothing_type can only be "shared" for "single" pathogen_type')
  }
  if (pathogen_str == 'single' & is.null(smoothing_obj$tau_priors$mean)){
    smoothing_obj$tau_priors$mean <- 0.0
    smoothing_obj$tau_priors$sd <- 1.0
    return(smoothing_obj)
  }

  # For independent structure, validate and adjust dimensions
  if (smoothing_obj$smoothing_type == "independent") {
    if (is.null(pathogen_names)) {
      stop("pathogen_names is required for 'independent' smoothing structure")
    }

    expected_dim <- length(pathogen_names)

    # If no priors were provided (priors_provided = 1), create dummy values with correct dimensions
    if (smoothing_obj$priors_provided == 1) {
      smoothing_obj$tau_priors$mean <- rep(0.0, expected_dim)
      smoothing_obj$tau_priors$sd <- rep(1.0, expected_dim)
    } else {
      # If priors were provided, validate and adjust dimensions
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

  # For shared structure with no priors, ensure values are 1-element arrays
  if (smoothing_obj$smoothing_type == "shared" && smoothing_obj$priors_provided == 1) {
    if (is.null(smoothing_obj$tau_priors$mean) || length(smoothing_obj$tau_priors$mean) == 0) {
      smoothing_obj$tau_priors$mean <- array(0.0, dim = 1)
      smoothing_obj$tau_priors$sd <- array(1.0, dim = 1)
    }
  }

  if (smoothing_obj$smoothing_type == "correlated" && smoothing_obj$priors_provided == 1) {
    if (is.null(smoothing_obj$tau_priors$mean) || length(smoothing_obj$tau_priors$mean) == 0) {
      smoothing_obj$tau_priors$mean <- numeric(0)
      smoothing_obj$tau_priors$sd <- numeric(0)
    }
  }

  return(smoothing_obj)
}

# Helper: Validate column exists in data
#' @noRd
check_column_exists <- function(data, col_name) {
  if (!col_name %in% colnames(data)) {
    stop(sprintf("Column '%s' not found in data", col_name))
  }
}

#' Helper: Detect if data is a time series class
#'
#' @param data input data object
#'
#' @returns logical indicating if data is a time series class
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
is_timeseries_class <- function(data) {
  ts_classes <- c("ts", "mts", "xts", "zoo", "zooreg", "tsibble", "tbl_ts",
                  "tibbletime", "tbl_time", "timeSeries", "irts", "tis")
  any(class(data) %in% ts_classes)
}

# Helper: Validate column is numeric
#' @srrstats {TS1.7} accommodate potential units data input
#' @noRd
check_column_numeric <- function(data, col_name) {
  if (!is.numeric(data[[col_name]]) && !inherits(data[[col_name]], "units")) {
    stop(sprintf("Column '%s' specified must be numeric or have units,
                 but is %s",
                 col_name, class(data[[col_name]])[1]))
  }
}

#' Check for missing data (NA values) in specified columns
#'
#' @param data dataframe to check
#' @param columns character vector of column names to check for missing data
#' @param context character string describing the context (e.g., "case_timeseries")
#'
#' @importFrom utils head
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G2.13, G2.14, G2.14a} check for missing data, error if found
#' @srrstats {G2.15, BS3.0} Data prep does not assume non-missingness
#'
#' @returns NULL (invisibly) if no missing data, otherwise throws an error
check_missing_data <- function(data, columns, context) {
  for (col in columns) {
    if (any(is.na(data[[col]]))) {
      na_count <- sum(is.na(data[[col]]))
      na_indices <- which(is.na(data[[col]]))

      # Show first few indices if there are many
      if (na_count > 5) {
        indices_msg <- paste0(
          "rows ",
          paste(utils::head(na_indices, 5), collapse = ", "),
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

#' Helper: Convert time series objects to tsibble
#'
#' @param ts_obj time series object (ts, mts, xts, zoo, etc.)
#'
#' @importFrom timetk tk_tbl
#'
#' @returns tsibble object in wide format
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
convert_ts_to_tsibble <- function(ts_obj) {

  # Handle xts and zoo objects using timetk (best practice for these classes)
  if (inherits(ts_obj, c("xts", "zoo", "zooreg"))) {

    # Check if timetk is available
    if (!requireNamespace("timetk", quietly = TRUE)) {
      cli::cli_abort(
        c(
          "The {.pkg timetk} package is required to handle
          {.cls {class(ts_obj)[1]}} objects.",
          "i" = "Install it with: install.packages('timetk')"
        )
      )
    }

    # Convert xts/zoo to tibble using timetk (preserves wide format)
    temp_tbl <- timetk::tk_tbl(ts_obj, preserve_index = TRUE, rename_index = "index")

    # Convert to tsibble
    temp_tsbl <- tryCatch({
      tsibble::as_tsibble(temp_tbl, index = .data$index)
    }, error = function(e) {
      cli::cli_abort(
        "Error converting {.cls {class(ts_obj)[1]}} to tsibble: {e$message}"
      )
    })

    return(temp_tsbl)
  }

  # Handle ts and mts objects using tsibble's built-in conversion
  if (inherits(ts_obj, c("ts", "mts"))) {
    temp_tsbl <- tryCatch({
      if (is.matrix(ts_obj)) {
        # Multivariate time series - keep wide format
        tsibble::as_tsibble(ts_obj, pivot_longer = FALSE)
      } else {
        # Univariate time series
        tsibble::as_tsibble(ts_obj)
      }
    }, error = function(e) {
      cli::cli_abort(
        "Error converting {.cls {class(ts_obj)[1]}} to tsibble: {e$message}"
      )
    })

    return(temp_tsbl)
  }

  # Handle tsibble objects (already in correct format)
  if (inherits(ts_obj, c("tsibble", "tbl_ts"))) {
    return(ts_obj)
  }

  # Try tsibble's as_tsibble for other time series classes
  temp_tsbl <- tryCatch({
    tsibble::as_tsibble(ts_obj)
  }, error = function(e) {
    # If that fails, provide a helpful error message
    cli::cli_abort(
      c(
        "Unable to convert time series object of class {.cls {class(ts_obj)}} to tsibble.",
        "i" = "Supported classes are: ts, mts, xts, zoo, zooreg, tsibble.",
        "i" = "For xts/zoo objects, the {.pkg timetk} package is required."
      )
    )
  })

  return(temp_tsbl)
}

#' Helper: Create and validate tsibble from various data formats
#'
#' @param data input dataset - can be data frame, tibble, data.table, or time series object
#' @param columns column names with data to include in the validated tsibble
#' @param time_col the column name with time data to validate (optional for ts classes)
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {TS1.1, TS1.2, TS1.3, TS1.4, TS1.5, TS1.6, TS2.0, TS2.1a} all
#'   three data ingest functions (in `pathogen_structure.R`) allow for any
#'   kind of timeseries class input or non-time-series input. This validation
#'   function converts them to a consistent format and run checks, including
#'   checking ordering (and reporting violations), and checking gaps in the
#'   time sequence (explicit missing values). NA values (implicit missing
#'   values) are not allowed due to the nature of the models.
create_validated_timeseries <- function(data, columns, time_col = NULL) {

  # Step 1: Handle time series objects
  if (is_timeseries_class(data)) {
    if (!is.null(time_col)) {
      cli::cli_warn(
        "The {.arg time} argument is ignored when data is a time series object.
        The time index will be extracted automatically."
      )
    }

    # Convert time series to tsibble, which will be in wide format
    temp_tsbl <- convert_ts_to_tsibble(data)

    # Get the time column name from the tsibble
    time_col <- as.character(tsibble::index_var(temp_tsbl))

    # Subset to required columns
    temp_tsbl <- temp_tsbl[, c(time_col, columns)]

    #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling
    #' @srrstats {TS1.7} accommodate units data input
    for (col in columns) {
      check_column_numeric(temp_tsbl, col)

      if (inherits(temp_tsbl[[col]], "units")) {
        temp_tsbl[[col]] <- as.numeric(temp_tsbl[[col]])
      }
    }

  } else {
    # Step 2: Handle data frame-like objects
    if (is.null(time_col)) {
      cli::cli_abort(
        "When {.arg time} is not specified, data must be a time series class
        object (ts, xts, zoo, tsibble, etc.)."
      )
    }

    # Convert to data frame if needed
    temp_df <- as.data.frame(data)

    # Subset to only the columns we need
    temp_df <- temp_df[, c(time_col, columns), drop = FALSE]

    #' @srrstats {G2.4, G2.4b} ensure consistent numeric handling
    #' @srrstats {TS1.7} accommodate units data input
    for (col in columns) {
      check_column_numeric(temp_df, col)

      if (inherits(temp_df[[col]], "units")) {
        temp_df[[col]] <- as.numeric(temp_df[[col]])
      }
    }

    check_missing_data(temp_df, columns, "input data")

    # Check ordering BEFORE creating tsibble (tsibble auto-sorts, hiding the issue)
    if (is.unsorted(temp_df[[time_col]])) {
      cli::cli_abort(
        "Time series is not ordered. Please sort the data in chronological order."
      )
    }

    # Create tsibble
    temp_tsbl <- tryCatch({
      tsibble::tsibble(temp_df, index = !!rlang::sym(time_col))
    }, error = function(e) {
      cli::cli_abort("Error creating tsibble: {e$message}")
    })
  }

  # Step 3: Validate tsibble properties (same for both paths)

  # Check for gaps
  if (tsibble::has_gaps(temp_tsbl)$.gaps) {
    cli::cli_abort(
      "Time series has gaps. All time points must be present with no missing
      time periods."
    )
  }

  # Check regularity
  if (!tsibble::is_regular(temp_tsbl)) {
    cli::cli_abort(
      "Time series is irregular. The model requires regularly spaced time
      intervals."
    )
  }

  return(temp_tsbl)
}

#' Validate generation interval distribution function
#'
#' @param gi_dist A function that takes a numeric vector and returns
#'   generation interval probabilities
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
validate_gi_dist <- function(gi_dist) {
  # Check if it's a function
  if (!is.function(gi_dist)) {
    cli::cli_abort(c(
      "{.arg gi_dist} must be a function",
      "x" = "You supplied a {.cls {class(gi_dist)}} object"
    ))
  }

  # Check if it accepts at least one argument
  if (length(formals(gi_dist)) < 1) {
    cli::cli_abort(c(
      "{.arg gi_dist} must accept at least one argument",
      "i" = "The function should accept {.arg x} as input"
    ))
  }

  # Test with a numeric vector
  test_input <- c(0, 1, 5, 10)
  test_output <- tryCatch(
    gi_dist(test_input),
    error = function(e) {
      cli::cli_abort(c(
        "{.arg gi_dist} failed when called with numeric input",
        "x" = "Error: {e$message}"
      ))
    }
  )

  # Check output is numeric
  if (!is.numeric(test_output)) {
    cli::cli_abort(c(
      "{.arg gi_dist} must return numeric values",
      "x" = "Function returned {.cls {class(test_output)}}"
    ))
  }

  # Check output has same length as input (vectorized)
  if (length(test_output) != length(test_input)) {
    cli::cli_abort(c(
      "{.arg gi_dist} must be vectorized (return same length as input)",
      "x" = "Input length: {length(test_input)}, Output length: {length(test_output)}"
    ))
  }

  # Check for non-negative values (probabilities should be >= 0)
  if (any(test_output < 0, na.rm = TRUE)) {
    cli::cli_abort(c(
      "{.arg gi_dist} must return non-negative values",
      "x" = "Found {sum(test_output < 0, na.rm = TRUE)} negative value{?s}"
    ))
  }

  # Check for NA/NaN/Inf
  if (any(!is.finite(test_output))) {
    cli::cli_warn(c(
      "{.arg gi_dist} returns non-finite values for some inputs",
      "i" = "Found {sum(!is.finite(test_output))} NA/NaN/Inf value{?s}"
    ))
  }

  invisible(NULL)
}

#' Validate pathogen combination arguments
#'
#' @param combination Vector of pathogen names or NULL
#' @param pathogen_names Vector of valid pathogen names from the model
#' @param arg_name Name of the argument being validated (for error messages).
#'   Must be either "numerator_combination" or "denominator_combination"
#' @noRd
validate_pathogen_combination <- function(combination, pathogen_names, arg_name) {

  if (is.null(combination)) {
    idx <- seq_along(pathogen_names)
    return(idx)
  }

  if (length(combination) == 0) {
    cli::cli_abort(c(
      "{.arg {arg_name}} cannot be empty",
      "i" = "Provide pathogen name{?s}, use {.code NULL} for default"
    ))
  }

  # Validate pathogen names
  invalid_names <- setdiff(combination, pathogen_names)
  if (length(invalid_names) > 0) {
    cli::cli_abort(c(
      "{.arg {arg_name}} contains invalid pathogen name{?s}",
      "x" = "Invalid: {.val {invalid_names}}",
      "i" = "Valid pathogen names: {.val {pathogen_names}}",
      NULL
    ))
  }

  if (!is.null(combination)) {
    idx <- match(combination, pathogen_names)
    return(idx)
  }
  invisible(NULL)
}

#' Validate R-hat Threshold Parameter
#'
#' Validates that rhat_threshold is a positive numeric value greater than or equal to 1
#'
#' @param rhat_threshold The R-hat convergence threshold
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_rhat_threshold <- function(rhat_threshold) {
  if (!is.numeric(rhat_threshold)) {
    cli::cli_abort("{.arg rhat_threshold} must be numeric, got {.cls {class(rhat_threshold)}}")
  }

  if (length(rhat_threshold) != 1) {
    cli::cli_abort("{.arg rhat_threshold} must be a single value, got length {length(rhat_threshold)}")
  }

  if (!is.finite(rhat_threshold)) {
    cli::cli_abort("{.arg rhat_threshold} must be a finite number")
  }

  if (rhat_threshold < 1.0) {
    cli::cli_abort(c(
      "{.arg rhat_threshold} must be greater than or equal to 1.0",
      "x" = "R-hat values are always >= 1 by definition",
      "i" = "Common thresholds: 1.01 (strict), 1.05 (moderate), 1.1 (lenient)"
    ))
  }

  # Warn if threshold seems unusual
  if (rhat_threshold > 2.0) {
    cli::cli_warn(c(
      "{.arg rhat_threshold} = {rhat_threshold} is unusually high",
      "i" = "This may not detect convergence issues effectively",
      "i" = "Consider using a value between 1.01 and 1.1"
    ))
  }

  invisible(NULL)
}

#' Validate Effective Sample Size Threshold Parameter
#'
#' Validates that eff_sample_threshold is a positive numeric value
#'
#' @param eff_sample_threshold The minimum effective sample size threshold
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_eff_sample_threshold <- function(eff_sample_threshold) {
  if (!is.numeric(eff_sample_threshold)) {
    cli::cli_abort("{.arg eff_sample_threshold} must be numeric, got {.cls {class(eff_sample_threshold)}}")
  }

  if (length(eff_sample_threshold) != 1) {
    cli::cli_abort("{.arg eff_sample_threshold} must be a single value, got length {length(eff_sample_threshold)}")
  }

  if (!is.finite(eff_sample_threshold)) {
    cli::cli_abort("{.arg eff_sample_threshold} must be a finite number")
  }

  if (eff_sample_threshold <= 0) {
    cli::cli_abort(c(
      "{.arg eff_sample_threshold} must be positive, got {eff_sample_threshold}",
      "i" = "Effective sample size represents the number of independent samples",
      "i" = "Common thresholds: 100 (minimum), 400 (good), 1000 (excellent)"
    ))
  }

  # Warn if threshold seems unusual
  if (eff_sample_threshold < 50) {
    cli::cli_warn(c(
      "{.arg eff_sample_threshold} = {eff_sample_threshold} is very low",
      "!" = "This may accept poorly sampled parameters",
      "i" = "Consider using a value of at least 100"
    ))
  }

  if (eff_sample_threshold > 5000) {
    cli::cli_warn(c(
      "{.arg eff_sample_threshold} = {eff_sample_threshold} is very high",
      "!" = "This threshold may be difficult to achieve",
      "i" = "Consider using a value between 100 and 1000"
    ))
  }

  invisible(NULL)
}
