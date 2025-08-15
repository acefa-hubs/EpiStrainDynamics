# =============================================================================
# SHARED HELPER FUNCTIONS
# =============================================================================

#' Calculate Summary Statistics for Posterior Samples
#'
#' Computes quantiles and proportion above threshold for a vector of posterior samples
#'
#' @param values Numeric vector of posterior samples
#' @param threshold Numeric threshold for calculating proportion above (default: 0)
#'
#' @importFrom stats quantile
#'
#' @return Data frame with columns: y (median), lb_50/ub_50 (50% CI), lb_95/ub_95 (95% CI), prop (proportion > threshold)
calc_stats <- function(values, threshold = 0) {
  quan <- stats::quantile(values, c(0.5, 0.025, 0.25, 0.75, 0.975))
  prop <- mean(values > threshold)

  data.frame(
    y = quan[[1]],
    lb_50 = quan[[3]],
    ub_50 = quan[[4]],
    lb_95 = quan[[2]],
    ub_95 = quan[[5]],
    prop = prop
  )
}

#' Extract Model Components from Fitted Model
#'
#' Extracts commonly used components from a fitted model object for analysis functions
#'
#' @param fitted_model Fitted model object containing fit, constructed_model, etc.
#' @return List containing extracted components: fit, pathogen_names, num_path,
#'   time_seq, time, num_days, days_per_knot, spline_degree, DOW, week_effect, dow_effect
get_model_components <- function(fitted_model) {
  list(
    fit = fitted_model$fit,
    pathogen_names = fitted_model$constructed_model$pathogen_names %||% NULL,
    num_path = length(fitted_model$constructed_model$pathogen_names %||% 1),
    time_seq = fitted_model$constructed_model$data$time_seq,
    time = fitted_model$constructed_model$data$time,
    num_days = length(fitted_model$constructed_model$data$time),
    days_per_knot = fitted_model$constructed_model$model_params$days_per_knot %||% NULL,
    spline_degree = fitted_model$constructed_model$model_params$spline_degree %||% NULL,
    DOW = fitted_model$constructed_model$model_params$DOW %||% NULL,
    week_effect = fitted_model$constructed_model$model_params$week_effect %||% NULL,
    dow_effect = fitted_model$constructed_model$dow_effect
  )
}

#' Expand Time Grid for Multiple Pathogens
#'
#' Expands a time grid to include all pathogen-time combinations (sequential
#' pairing, not Cartesian product)
#'
#' @param time_grid Data frame from with time_idx columns
#' @param pathogen_names Character vector of pathogen names
#' @return Data frame with columns: pathogen, pathogen_idx, time_idx
expand_pathogen_grid <- function(time_grid, pathogen_names) {
  num_days <- nrow(time_grid)
  num_path <- length(pathogen_names)

  data.frame(
    pathogen = rep(pathogen_names, each = num_days),
    pathogen_idx = rep(seq_along(pathogen_names), each = num_days),
    time_idx = rep(time_grid$time_idx, times = num_path)
  )
}

#' Predict B-spline Basis Matrix
#'
#' Creates B-spline basis matrix for transforming spline coefficients to time series
#'
#' @param time_seq Numeric vector of time points
#' @param days_per_knot Integer number of days between knots
#' @param spline_degree Integer degree of B-splines
#' @return Matrix B_true for transforming spline coefficients
#' @importFrom splines bs
#' @importFrom stats predict
predict_B_true <- function(time_seq, days_per_knot, spline_degree) {
  X <- as.numeric(time_seq)
  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)

  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1),
                        knots = knots[2:(length(knots)-1)],
                        degree = spline_degree,
                        intercept = TRUE)
  t(stats::predict(B_true, X))
}

#' Transform Posterior Samples Using Splines (Single Pathogen)
#'
#' Transforms spline coefficients to time series for single pathogen models
#'
#' @param post Posterior samples object containing spline coefficients
#' @param B_true B-spline basis matrix from predict_B_true()
#' @param num_days Integer number of time points
#' @return Array of transformed posterior samples [samples, time]
transform_posterior_single <- function(post, B_true, num_days) {
  a <- array(data = NA, dim = c(nrow(post$a), num_days))
  for(k in 1:nrow(post$a)) {
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }
  a
}

#' Transform Posterior Samples Using Splines (Multiple Pathogens)
#'
#' Transforms spline coefficients to time series for multiple pathogen models
#'
#' @param post Posterior samples object containing spline coefficients
#' @param B_true B-spline basis matrix from predict_B_true()
#' @param num_path Integer number of pathogens
#' @param num_days Integer number of time points
#' @return Array of transformed posterior samples [samples, pathogens, time]
transform_posterior_multi <- function(post, B_true, num_path, num_days) {
  a <- array(data = NA, dim = c(nrow(post$a), num_path, num_days))
  for(j in 1:num_path) {
    for(k in 1:nrow(post$a)) {
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }
  a
}

# =============================================================================
# COMPUTATION ENGINE
# =============================================================================

#' Generic Computation Engine for Multiple Pathogen Analysis
#'
#' High-level function that coordinates analysis for multiple pathogen models
#'
#' @param fitted_model Fitted model object
#' @param start_idx Integer starting time index for analysis
#' @param measure Character string specifying metric ("incidence", "growth_rate", "Rt")
#' @param threshold Numeric threshold for proportion calculations (default: 0)
#' @param use_splines Logical indicating whether to use spline transformation
#' @param ... Additional arguments passed to calculation functions
#'  (e.g., tau_max, gi_dist for Rt, and dow for incidence)
#' @return Data frame with results for individual pathogens and totals
#' @importFrom rstan extract
#' @importFrom dplyr bind_rows arrange
#' @examples
#' \dontrun{
#' results <- compute_multi_pathogen(fitted_model, 1, "incidence")
#' }
compute_multi_pathogen <- function(fitted_model, start_idx, measure,
                                   threshold = 0, use_splines = FALSE,
                                   ...) {

  components <- get_model_components(fitted_model)
  post <- rstan::extract(components$fit)

  # Transform data if using splines
  if (use_splines) {
    B_true <- predict_B_true(components$time_seq, components$days_per_knot,
                             components$spline_degree)
    a <- transform_posterior_multi(post, B_true, components$num_path,
                                   components$num_days)
  } else {
    a <- post$a
  }

  selection_index <- start_idx:components$num_days
  time_grid <- data.frame(time_idx = selection_index)
  extra_args <- list(...)

  # Individual pathogen results
  pathogen_grid <- expand_pathogen_grid(time_grid, components$pathogen_names)

  calc_individual_pathogen_fn <- switch(
    measure,
    "incidence" = calc_incidence_individual,
    "growth_rate" = calc_growth_individual,
    "Rt" = calc_rt_individual,
    stop("Unknown metric: ", measure)
  )
  pathogen_results <- calc_wrapper(pathogen_grid, pathogen_grid$time_idx,
                                   pathogen_grid$pathogen_idx,
                                   calc_individual_pathogen_fn,
                                   a, post, components, extra_args, threshold)

  calc_total_pathogens_fn <- switch(
    measure,
    "incidence" = calc_incidence_total,
    "growth_rate" = calc_growth_total,
    "Rt" = calc_rt_total,
    stop("Unknown metric: ", measure)
  )
  total_results <- calc_wrapper(time_grid, time_grid$time_idx,
                                pathogen_idx_col = rep(list(NULL), nrow(time_grid)),
                                calc_total_pathogens_fn,
                                a, post, components, extra_args, threshold)
  total_results$pathogen <- "Total"

  measure_out <- dplyr::bind_rows(pathogen_results, total_results) |>
    dplyr::arrange(pathogen != "Total", pathogen) |>
    cbind(time = components$time[selection_index])

  out <- list(measure = measure_out,
              fit = fitted_model$fit,
              constructed_model = fitted_model$constructed_model)
}

#' Generic Computation Engine for Single Pathogen Analysis
#'
#' High-level function that coordinates analysis for single pathogen models
#'
#' @param fitted_model Fitted model object
#' @param start_idx Integer starting time index for analysis
#' @param measure Character string specifying metric ("incidence", "growth_rate", "Rt")
#' @param threshold Numeric threshold for proportion calculations (default: 0)
#' @param use_splines Logical indicating whether to use spline transformation
#' @param ... Additional arguments passed to calculation functions (e.g., tau_max, gi_dist for Rt)
#' @return Data frame with analysis results
#' @importFrom rstan extract
#' @examples
#' \dontrun{
#' results <- compute_single_pathogen(fitted_model, 1, "growth_rate")
#' }
compute_single_pathogen <- function(fitted_model, start_idx, measure,
                                    threshold = 0, use_splines = FALSE,
                                    ...) {

  components <- get_model_components(fitted_model)
  post <- rstan::extract(components$fit)

  # Transform data if using splines
  if (use_splines) {
    B_true <- predict_B_true(components$time_seq, components$days_per_knot, components$spline_degree)
    a <- transform_posterior_single(post, B_true, components$num_days)
  } else {
    a <- post$a
  }

  extra_args <- list(...)

  # Create results
  selection_index <- start_idx:components$num_days
  time_grid <- data.frame(time_idx = selection_index)

  calc_single_pathogen_fn <- switch(
    measure,
    "incidence" = calc_incidence_single,
    "growth_rate" = calc_growth_single,
    "Rt" = calc_rt_single,
    stop("Unknown metric: ", measure)
  )
  results <- calc_wrapper(time_grid, time_grid$time_idx,
                          pathogen_idx_col = rep(list(NULL), nrow(time_grid)),
                          calc_single_pathogen_fn,
                          a, post, components, extra_args, threshold)

  measure <- cbind(results,
                   time = components$time[selection_index])

  out <- list(measure = measure,
              fit = fitted_model$fit,
              constructed_model = fitted_model$constructed_model)
}

#' Proportion analysis for multiple pathogens
#'
#' High-level function to coordinates proportion analysis for multi pathogen models
#'
#' @param fitted_model Fitted model object
#' @param numerator_idx Integer with index of pathogens to place in numerator
#' @param threshold Numeric threshold for proportion calculations (default: 0)
#' @param use_splines Logical indicating whether to use spline transformation
#' @param ... Additional arguments passed to calculation functions (e.g., tau_max, gi_dist for Rt)
#' @return Data frame with analysis results
#' @importFrom rstan extract
#' @examples
#' \dontrun{
#' results <- compute_single_pathogen(fitted_model, 1, "growth_rate")
#' }
compute_proportions <- function(fitted_model, numerator_idx,
                                threshold = 0, use_splines = FALSE,
                                ...) {

  components <- get_model_components(fitted_model)
  post <- rstan::extract(components$fit)

  # Transform data if using splines
  if (use_splines) {
    B_true <- predict_B_true(components$time_seq, components$days_per_knot, components$spline_degree)
    a <- transform_posterior_multi(post, B_true, components$num_path, components$num_days)
  } else {
    a <- post$a
  }

  extra_args <- list(...)

  # Create results
  time_grid <- data.frame(time_idx = 1:components$num_days)

  results <- calc_wrapper(time_grid, time_grid$time_idx,
                          pathogen_idx_col = rep(list(numerator_idx), nrow(time_grid)),
                          calc_proportion,
                          a, post, components, extra_args, threshold)
  measure <- cbind(results, components$time)

  out <- list(measure = measure,
              fit = fitted_model$fit,
              constructed_model = fitted_model$constructed_model)

  return(out)
}

#' Unified Calculation Wrapper
#'
#' Wrapper function that applies calculation functions to data and computes summary statistics
#'
#' @param df Data frame containing time grid or pathogen grid
#' @param time_idx_col Vector of time indices
#' @param pathogen_idx_col Vector of pathogen indices (or list of NULLs)
#' @param calc_fn Calculation function to apply
#' @param a Array of posterior samples
#' @param post Posterior samples object
#' @param components Model components from get_model_components()
#' @param extra_args List of additional arguments for calc_fn
#' @param threshold Numeric threshold for summary statistics
#' @return Data frame with expanded summary statistics
#' @importFrom purrr map2
calc_wrapper <- function (df, time_idx_col, pathogen_idx_col, calc_fn,
                          a, post, components, extra_args, threshold) {

  # Calculate statistics for each time/pathogen combination
  stats_list <- purrr::map2(time_idx_col, pathogen_idx_col,
                            function(t_idx, p_idx) {
                              values <- do.call(
                                calc_fn, c(list(a = a,
                                                time_idx = t_idx,
                                                pathogen_idx = p_idx,
                                                post = post,
                                                components = components),
                                           extra_args))
                              calc_stats(values, threshold = threshold)
                            })

  # Convert list of data frames to single data frame
  stats_df <- do.call(rbind, stats_list)

  # Combine with original data frame (excluding time_idx column)
  result_df <- cbind(df[, !names(df) %in% "time_idx", drop = FALSE], stats_df)

  return(result_df)
}
