#' Validate MCMC Chain Parameter
#'
#' Validates that n_chain is a positive integer
#'
#' @param n_chain The number of MCMC chains
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_n_chain <- function(n_chain) {
  if (!is.numeric(n_chain)) {
    cli::cli_abort("{.arg n_chain} must be numeric, got {.cls {class(n_chain)}}")
  }

  if (length(n_chain) != 1) {
    cli::cli_abort("{.arg n_chain} must be a single value, got length {length(n_chain)}")
  }

  if (!is.finite(n_chain)) {
    cli::cli_abort("{.arg n_chain} must be a finite number")
  }

  if (n_chain != as.integer(n_chain)) {
    cli::cli_abort("{.arg n_chain} must be a whole number, got {n_chain}")
  }

  if (n_chain <= 0) {
    cli::cli_abort("{.arg n_chain} must be positive, got {n_chain}")
  }

  invisible(NULL)
}

#' Validate MCMC Iteration Parameter
#'
#' Validates that n_iter is a positive integer
#'
#' @param n_iter The number of MCMC iterations
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_n_iter <- function(n_iter) {
  if (!is.numeric(n_iter)) {
    cli::cli_abort("{.arg n_iter} must be numeric, got {.cls {class(n_iter)}}")
  }

  if (length(n_iter) != 1) {
    cli::cli_abort("{.arg n_iter} must be a single value, got length {length(n_iter)}")
  }

  if (!is.finite(n_iter)) {
    cli::cli_abort("{.arg n_iter} must be a finite number")
  }

  if (n_iter != as.integer(n_iter)) {
    cli::cli_abort("{.arg n_iter} must be a whole number, got {n_iter}")
  }

  if (n_iter <= 0) {
    cli::cli_abort("{.arg n_iter} must be positive, got {n_iter}")
  }

  invisible(NULL)
}

#' Validate MCMC Warmup Parameter
#'
#' Validates that n_warmup is a non-negative integer
#'
#' @param n_warmup The number of warmup iterations
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_n_warmup <- function(n_warmup) {
  if (!is.numeric(n_warmup)) {
    cli::cli_abort("{.arg n_warmup} must be numeric, got {.cls {class(n_warmup)}}")
  }

  if (length(n_warmup) != 1) {
    cli::cli_abort("{.arg n_warmup} must be a single value, got length {length(n_warmup)}")
  }

  if (!is.finite(n_warmup)) {
    cli::cli_abort("{.arg n_warmup} must be a finite number")
  }

  if (n_warmup != as.integer(n_warmup)) {
    cli::cli_abort("{.arg n_warmup} must be a whole number, got {n_warmup}")
  }

  if (n_warmup < 0) {
    cli::cli_abort("{.arg n_warmup} cannot be negative, got {n_warmup}")
  }

  invisible(NULL)
}

#' Validate MCMC Thinning Parameter
#'
#' Validates that thin is a positive integer
#'
#' @param thin The thinning interval
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_thin <- function(thin) {
  if (!is.numeric(thin)) {
    cli::cli_abort("{.arg thin} must be numeric, got {.cls {class(thin)}}")
  }

  if (length(thin) != 1) {
    cli::cli_abort("{.arg thin} must be a single value, got length {length(thin)}")
  }

  if (!is.finite(thin)) {
    cli::cli_abort("{.arg thin} must be a finite number")
  }

  if (thin != as.integer(thin)) {
    cli::cli_abort("{.arg thin} must be a whole number, got {thin}")
  }

  if (thin <= 0) {
    cli::cli_abort("{.arg thin} must be positive, got {thin}")
  }

  invisible(NULL)
}

#' Validate adapt_delta Parameter
#'
#' Validates that adapt_delta is between 0 and 1
#'
#' @param adapt_delta The target acceptance probability
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_adapt_delta <- function(adapt_delta) {
  if (!is.numeric(adapt_delta)) {
    cli::cli_abort("{.arg adapt_delta} must be numeric, got {.cls {class(adapt_delta)}}")
  }

  if (length(adapt_delta) != 1) {
    cli::cli_abort("{.arg adapt_delta} must be a single value, got length {length(adapt_delta)}")
  }

  if (!is.finite(adapt_delta)) {
    cli::cli_abort("{.arg adapt_delta} must be a finite number")
  }

  if (adapt_delta <= 0 || adapt_delta >= 1) {
    cli::cli_abort("{.arg adapt_delta} must be between 0 and 1, got {adapt_delta}")
  }

  invisible(NULL)
}

#' Validate seed Parameter
#'
#' Validates that seed is NULL or a positive integer
#'
#' @param seed The random seed
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }

  if (!is.numeric(seed)) {
    cli::cli_abort("{.arg seed} must be numeric or NULL, got {.cls {class(seed)}}")
  }

  if (length(seed) != 1) {
    cli::cli_abort("{.arg seed} must be a single value, got length {length(seed)}")
  }

  if (!is.finite(seed)) {
    cli::cli_abort("{.arg seed} must be a finite number")
  }

  if (seed != as.integer(seed)) {
    cli::cli_abort("{.arg seed} must be a whole number, got {seed}")
  }

  if (seed <= 0) {
    cli::cli_abort("{.arg seed} must be positive, got {seed}")
  }

  invisible(NULL)
}

#' Validate verbose Parameter
#'
#' Validates that verbose is a logical value
#'
#' @param verbose The verbosity flag
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_verbose <- function(verbose) {
  if (!is.logical(verbose)) {
    cli::cli_abort("{.arg verbose} must be logical (TRUE/FALSE), got {.cls {class(verbose)}}")
  }

  if (length(verbose) != 1) {
    cli::cli_abort("{.arg verbose} must be a single value, got length {length(verbose)}")
  }

  if (is.na(verbose)) {
    cli::cli_abort("{.arg verbose} cannot be NA")
  }

  invisible(NULL)
}

#' Validate multi_cores Parameter
#'
#' Validates that multi_cores is a logical value
#'
#' @param multi_cores The parallel processing flag
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_multi_cores <- function(multi_cores) {
  if (!is.logical(multi_cores)) {
    cli::cli_abort("{.arg multi_cores} must be logical (TRUE/FALSE), got {.cls {class(multi_cores)}}")
  }

  if (length(multi_cores) != 1) {
    cli::cli_abort("{.arg multi_cores} must be a single value, got length {length(multi_cores)}")
  }

  if (is.na(multi_cores)) {
    cli::cli_abort("{.arg multi_cores} cannot be NA")
  }

  invisible(NULL)
}
