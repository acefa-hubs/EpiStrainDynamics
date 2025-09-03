#' Case-insensitive argument matching
#'
#' A helper function to perform case-insensitive argument matching, similar to
#' `match.arg()` but ignoring case differences. This function converts both the
#' argument and choices to lowercase before matching.
#'
#' @param arg The argument to match. Can be a character vector or missing.
#' @param choices Character vector of valid choices (should be in lowercase).
#' @param several.ok Logical; if TRUE, multiple matches are allowed.
#' @param arg_name Optional character string giving the name of the argument
#'   for better error messages.
#'
#' @return The matched argument (in lowercase).
match_arg_case_insensitive <- function(arg, choices, several.ok = FALSE, arg_name = NULL) {

  # Handle missing argument case (use first choice as default)
  if (missing(arg)) {
    return(choices[1L])
  }

  # Handle case where arg is the full choices vector (default parameter case)
  if (length(arg) > 1 && identical(sort(tolower(arg)), sort(tolower(choices)))) {
    return(choices[1L])
  }

  # Convert to lowercase for matching
  arg_lower <- tolower(arg)
  choices_lower <- tolower(choices)

  # Perform the matching
  matched_indices <- pmatch(arg_lower, choices_lower)

  # Handle multiple matches
  if (length(matched_indices) > 1L && !several.ok) {
    arg_display <- if (is.null(arg_name)) "argument" else paste0("'", arg_name, "'")
    stop(paste0(arg_display, " must be of length 1"), call. = FALSE)
  }

  # Handle no matches
  if (any(is.na(matched_indices))) {
    arg_display <- if (is.null(arg_name)) "argument" else paste0("'", arg_name, "'")
    invalid_args <- arg[is.na(matched_indices)]
    stop(paste0(arg_display, " should be one of: ",
                paste(paste0("'", choices, "'"), collapse = ", "),
                ". Got: ", paste(paste0("'", invalid_args, "'"), collapse = ", ")),
         call. = FALSE)
  }

  # Return the matched choice(s) in original case
  return(choices[matched_indices])
}

#' Apply case-insensitive matching to multiple arguments
#'
#' A convenience wrapper around `match_arg_case_insensitive()` for handling
#' multiple arguments at once. Useful when you have several string parameters
#' that all need case-insensitive matching.
#'
#' @param arg_list Named list where names are argument names and values are
#'   the arguments to match
#' @param choices_list Named list where names correspond to argument names
#'   and values are character vectors of valid choices
#' @param several.ok Logical; if TRUE, multiple matches are allowed for all args
#'
#' @return Named list with matched arguments (in lowercase)
match_args_case_insensitive <- function(arg_list, choices_list, several.ok = FALSE) {

  # Validate inputs
  if (!is.list(arg_list) || !is.list(choices_list)) {
    stop("Both arg_list and choices_list must be lists", call. = FALSE)
  }

  if (!all(names(arg_list) %in% names(choices_list))) {
    missing_choices <- setdiff(names(arg_list), names(choices_list))
    stop(paste0("Missing choices for arguments: ",
                paste(paste0("'", missing_choices, "'"), collapse = ", ")),
         call. = FALSE)
  }

  # Apply case-insensitive matching to each argument
  result <- lapply(names(arg_list), function(arg_name) {
    match_arg_case_insensitive(
      arg = arg_list[[arg_name]],
      choices = choices_list[[arg_name]],
      several.ok = several.ok,
      arg_name = arg_name
    )
  })

  names(result) <- names(arg_list)
  return(result)
}

#' Validate that input is a positive whole number
#'
#' @param value The value to validate
#' @param arg_name Name of the argument being validated (for error messages)
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
#'
validate_positive_whole_number <- function(value, arg_name) {

  # Check that value is numeric
  if (!is.numeric(value)) {
    stop(paste("Argument", paste0("'", arg_name, "'"), "must be numeric"),
         call. = FALSE)
  }

  # Check that value is a single value (length = 1)
  if (length(value) != 1) {
    stop(paste("Argument", paste0("'", arg_name, "'"), "must be a single value"),
         call. = FALSE)
  }

  # Check that value is finite (not NA, NULL, Inf, -Inf, NaN)
  if (!is.finite(value)) {
    stop(paste("Argument", paste0("'", arg_name, "'"), "must be a finite number"),
         call. = FALSE)
  }

  # Check that value is a whole number
  if (value != as.integer(value)) {
    stop(paste("Argument", paste0("'", arg_name, "'"), "must be a whole number"),
         call. = FALSE)
  }

  # Check that value is positive
  if (value <= 0) {
    stop(paste("Argument", paste0("'", arg_name, "'"), "must be a positive number"),
         call. = FALSE)
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
#' @return Invisible NULL if validation passes, otherwise throws an error
validate_class_inherits <- function(obj, class_names, require_all = TRUE) {
  # Ensure class_names is a character vector
  if (!is.character(class_names) || length(class_names) == 0) {
    stop("class_names must be a non-empty character vector", call. = FALSE)
  }

  # Check inheritance for each class
  inheritance_results <- sapply(class_names, function(cls) inherits(obj, cls))

  if (require_all) {
    # Must inherit from ALL classes
    if (!all(inheritance_results)) {
      missing_classes <- class_names[!inheritance_results]
      if (length(class_names) == 1) {
        stop(paste("Input must be of class", class_names,
                   "but got class:", paste(class(obj), collapse = ", ")),
             call. = FALSE)
      } else {
        stop(paste("Input must inherit from all classes:",
                   paste(class_names, collapse = ", "),
                   "\nMissing classes:", paste(missing_classes, collapse = ", "),
                   "\nActual classes:", paste(class(obj), collapse = ", ")),
             call. = FALSE)
      }
    }
  } else {
    # Must inherit from at least ONE class
    if (!any(inheritance_results)) {
      stop(paste("Input must inherit from at least one of:",
                 paste(class_names, collapse = ", "),
                 "\nActual classes:", paste(class(obj), collapse = ", ")),
           call. = FALSE)
    }
  }

  invisible(NULL)
}

#' Validate input vectors have matching lengths
#'
#' @param ... Named vectors to check for matching lengths
#' @param .error_call The calling function for better error messages
#'
validate_matching_lengths <- function(..., .error_call = rlang::caller_env()) {
  vectors <- list(...)
  vector_names <- names(vectors)
  lengths <- lengths(vectors)

  if (length(unique(lengths)) > 1) {
    length_info <- paste(paste(vector_names, lengths, sep = " (length "), ")", collapse = ", ")
    rlang::abort(
      message = paste("All input vectors must have the same length.",
                      "Found:", length_info),
      call = .error_call
    )
  }
}

#' Validate that a list of vectors all have the same length
#'
#' @param vector_list Named list of vectors
#' @param list_name Name of the list parameter for error messages
#' @param reference_length Expected length (optional)
#' @param .error_call The calling function for better error messages
#'
validate_list_vector_lengths <- function(vector_list, list_name, reference_length = NULL, .error_call = rlang::caller_env()) {
  if (length(vector_list) == 0) {
    rlang::abort(
      message = paste(list_name, "cannot be empty"),
      call = .error_call
    )
  }

  lengths <- lengths(vector_list)
  vector_names <- names(vector_list)

  if (length(unique(lengths)) > 1) {
    length_info <- paste(paste(vector_names, lengths, sep = " (length "), ")", collapse = ", ")
    rlang::abort(
      message = paste("All vectors in", list_name, "must have the same length.",
                      "Found:", length_info),
      call = .error_call
    )
  }

  if (!is.null(reference_length) && lengths[1] != reference_length) {
    rlang::abort(
      message = paste("Vectors in", list_name, "must have length", reference_length,
                      "but found length", lengths[1]),
      call = .error_call
    )
  }
}
