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
    cli::cli_abort("{.arg multi_cores} must be logical (TRUE/FALSE),
                   got {.cls {class(multi_cores)}}")
  }

  if (length(multi_cores) != 1) {
    cli::cli_abort("{.arg multi_cores} must be a single value, got length
                   {length(multi_cores)}")
  }

  if (is.na(multi_cores)) {
    cli::cli_abort("{.arg multi_cores} cannot be NA")
  }

  invisible(NULL)
}

#' Validate MCMC Parameters Collectively
#'
#' Performs cross-parameter validation and checks for extreme values that may
#' lead to long computation times or resource issues.
#'
#' @param n_iter The number of MCMC iterations
#' @param n_warmup The number of warmup iterations
#' @param n_chain The number of MCMC chains
#' @param thin The thinning interval
#' @param seed Random seed (optional)
#' @param suppress_warnings Logical; if TRUE, suppresses advisory warnings
#'
#' @noRd
#' @srrstats {BS2.6} Check that values for computational parameters lie within
#'   plausible ranges, including cross-parameter validation and warnings for
#'   extreme values that may cause long computation times
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
#'   or issues warnings
validate_mcmc_params_collective <- function(n_iter, n_warmup,
                                            n_chain, thin = 1, seed = NULL,
                                            suppress_warnings = FALSE) {
  # BS2.6: Cross-parameter validation - warmup must be less than total iterations
  if (n_warmup >= n_iter) {
    cli::cli_abort(
      c(
        "{.arg n_warmup} must be less than {.arg n_iter}",
        "x" = "Got {.arg n_warmup} = {n_warmup} and {.arg n_iter} = {n_iter}"
      )
    )
  }

  if (!suppress_warnings) {
    # BS2.6: Validate that thinning doesn't eliminate too many samples
    effective_samples <- floor((n_iter - n_warmup) / thin) * n_chain
    if (effective_samples < 100) {
      cli::cli_warn(
        c(
          "Very few effective samples will be retained after thinning",
          "i" = "With {.arg n_iter} = {n_iter}, {.arg n_warmup} = {n_warmup}, {.arg thin} = {thin}, and {.arg n_chain} = {n_chain}",
          "i" = "Only {effective_samples} samples will be retained",
          "i" = "Consider reducing {.arg thin} or increasing {.arg n_iter}"
        )
      )
    }

    # BS2.6: Warn about extremely large iteration counts (long computation time)
    if (n_iter > 20000) {
      cli::cli_warn(
        c(
          "Large {.arg n_iter} may result in long computation times",
          "i" = "Requested {n_iter} iterations",
          "i" = "Consider starting with fewer iterations to assess computational burden"
        )
      )
    }

    # BS2.6: Warn if requesting more chains than available cores
    n_cores <- parallel::detectCores(logical = FALSE)
    if (!is.na(n_cores) && n_chain > n_cores) {
      cli::cli_warn(
        c(
          "{.arg n_chain} exceeds available CPU cores",
          "i" = "Requested {n_chain} chains but only {n_cores} cores detected",
          "i" = "This may slow computation. Consider reducing {.arg n_chain} to {n_cores} or fewer"
        )
      )
    }
  }

  # BS2.6: Validate seed if provided — errors and coercion warnings kept
  # regardless of suppress_warnings since these indicate data issues
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      cli::cli_abort(
        c(
          "{.arg seed} must be a single numeric value or NULL",
          "x" = "Got {.cls {class(seed)}} of length {length(seed)}"
        )
      )
    }
    if (!is.finite(seed)) {
      cli::cli_abort("{.arg seed} must be a finite number")
    }
    if (seed != as.integer(seed)) {
      cli::cli_warn(
        c(
          "{.arg seed} is not an integer and will be coerced",
          "i" = "Got {seed}, will use {as.integer(seed)}"
        )
      )
    }
  }

  invisible(NULL)
}

#' Validate suppress_warnings Parameter
#'
#' Validates that suppress_warnings is a single logical value
#'
#' @param suppress_warnings The suppress_warnings parameter value
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#' @srrstats {BS2.14} validates suppress_warnings parameter
#'
#' @return NULL (invisibly) if validation passes, otherwise stops with error
validate_suppress_warnings <- function(suppress_warnings) {
  if (!is.logical(suppress_warnings)) {
    cli::cli_abort("{.arg suppress_warnings} must be logical (TRUE or FALSE),
                   got {.cls {class(suppress_warnings)}}")
  }

  if (length(suppress_warnings) != 1) {
    cli::cli_abort("{.arg suppress_warnings} must be a single value, got
                   length {length(suppress_warnings)}")
  }

  if (is.na(suppress_warnings)) {
    cli::cli_abort("{.arg suppress_warnings} cannot be NA")
  }

  invisible(NULL)
}
