#' Validate class of argument is correct
#'
#' @param obj argument to check
#' @param class_name class name the object should be
#'
validate_class_inherits <- function(obj, class_name) {
  if (!inherits(obj, class_name)) {
    stop(paste("Input must be of class", class_name,
               "but got class:", class(obj)[1]),
         call. = FALSE)
  }
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
