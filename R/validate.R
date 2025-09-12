#' Diagnose model convergence and fit
#'
#' @param fitted_model fitted model object
#' @param rhat_threshold R-hat threshold for convergence (default 1.1)
#' @param eff_sample_threshold effective sample size threshold
#'
#' @return list of diagnostic information
#' @export
diagnose_model <- function(fitted_model, rhat_threshold = 1.1, eff_sample_threshold = 100) {

  if (!inherits(fitted_model, "EpiStrainDynamics.fit")) {
    stop("fitted_model must be an EpiStrainDynamics.fit object")
  }

  # Get model summary
  fit_summary <- rstan::summary(fitted_model$fit)$summary

  # Check R-hat values
  rhat_values <- fit_summary[, "Rhat"]
  rhat_issues <- names(rhat_values)[rhat_values > rhat_threshold & !is.na(rhat_values)]

  # Check effective sample sizes
  n_eff_values <- fit_summary[, "n_eff"]
  eff_sample_issues <- names(n_eff_values)[n_eff_values < eff_sample_threshold & !is.na(n_eff_values)]

  # Overall convergence assessment
  convergence <- length(rhat_issues) == 0 && length(eff_sample_issues) == 0

  # Additional diagnostics
  max_rhat <- max(rhat_values, na.rm = TRUE)
  min_neff <- min(n_eff_values, na.rm = TRUE)

  diagnostics <- list(
    convergence = convergence,
    rhat_issues = rhat_issues,
    eff_sample_issues = eff_sample_issues,
    max_rhat = max_rhat,
    min_neff = min_neff,
    summary = fit_summary
  )

  # Print summary
  cat("Model Convergence Diagnostics\n")
  cat("=============================\n")
  cat("Overall convergence:", ifelse(convergence, "GOOD", "ISSUES DETECTED"), "\n")
  cat("Maximum R-hat:", round(max_rhat, 3), "\n")
  cat("Minimum n_eff:", round(min_neff, 0), "\n")

  return(invisible(diagnostics))
}

# Example usage:
# diagnostics <- diagnose_model(fitted_model)
# if (!diagnostics$convergence) {
#   warning("Model has convergence issues")
# }


# Main validation function users can call
validate_fit <- function(model_fit, method = c("convergence", "ppc", "loo", "all")) {
  method <- match.arg(method)

  results <- list()

  if (method %in% c("convergence", "all")) {
    results$convergence <- check_convergence(model_fit)
  }

  if (method %in% c("ppc", "all")) {
    results$posterior_predictive <- posterior_predictive_check(model_fit, model_fit$data)
  }

  if (method %in% c("loo", "all")) {
    results$cross_validation <- validate_model_loo(model_fit)
  }

  # Overall assessment
  results$validation_summary <- summarize_validation(results)

  class(results) <- "validation_results"
  return(results)
}

# Print method for validation results
#' @export
print.validation_results <- function(x) {
  cat("Model Validation Results\n")
  cat("========================\n\n")

  if (!is.null(x$convergence)) {
    cat("Convergence: ", ifelse(x$convergence$converged, "✓ PASSED", "✗ FAILED"), "\n")
  }

  if (!is.null(x$posterior_predictive)) {
    cat("Posterior Predictive Check: p-value =", round(x$posterior_predictive$p_value, 3), "\n")
  }

  if (!is.null(x$cross_validation)) {
    cat("Cross-validation: LOOIC =", round(x$cross_validation$looic, 2), "\n")
  }
}








#' Extract smoothing parameters from P-spline models
#'
#' @param fitted_model fitted P-spline model object
#'
#' @return list of smoothing parameter estimates
#' @export
extract_smoothing_parameters <- function(fitted_model) {

  if (!inherits(fitted_model, "EpiStrainDynamics.fit")) {
    stop("fitted_model must be an EpiStrainDynamics.fit object")
  }

  # Check if model is P-spline
  if (!inherits(fitted_model, "ps") && !inherits(fitted_model, "ps_single")) {
    stop("extract_smoothing_parameters only applies to P-spline models")
  }

  # Extract smoothing parameters
  posterior <- extract_posterior(fitted_model, pars = c("tau", "sigma_spline"))

  # Calculate summary statistics
  results <- list()

  if ("tau" %in% names(posterior)) {
    tau_samples <- posterior$tau
    results$tau <- list(
      mean = mean(tau_samples),
      median = median(tau_samples),
      sd = sd(tau_samples),
      quantiles = quantile(tau_samples, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      samples = tau_samples
    )
  }

  if ("sigma_spline" %in% names(posterior)) {
    sigma_samples <- posterior$sigma_spline
    results$sigma_spline <- list(
      mean = mean(sigma_samples),
      median = median(sigma_samples),
      sd = sd(sigma_samples),
      quantiles = quantile(sigma_samples, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
      samples = sigma_samples
    )
  }

  return(results)
}

# Example usage:
# smoothing_params <- extract_smoothing_parameters(ps_fitted_model)
# hist(smoothing_params$tau$samples, main = "Posterior Distribution of Tau")
