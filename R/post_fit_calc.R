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
#' @return List containing extracted components: fit, pathogen_names, num_path, time,
#'   time_labels, num_days, days_per_knot, spline_degree, DOW, week_effect, dow_effect
extract_model_components <- function(fitted_model) {
  list(
    fit = fitted_model$fit,
    pathogen_names = fitted_model$constructed_model$pathogen_names %||% NULL,
    num_path = length(fitted_model$constructed_model$pathogen_names %||% 1),
    time = fitted_model$constructed_model$data$time,
    time_labels = fitted_model$constructed_model$data$time,
    num_days = length(fitted_model$constructed_model$data$time),
    days_per_knot = fitted_model$constructed_model$model_params$days_per_knot %||% NULL,
    spline_degree = fitted_model$constructed_model$model_params$spline_degree %||% NULL,
    DOW = fitted_model$constructed_model$model_params$DOW %||% NULL,
    week_effect = fitted_model$constructed_model$model_params$week_effect %||% NULL,
    dow_effect = fitted_model$constructed_model$dow_effect
  )
}

#' Create Time Grid for Analysis
#'
#' Creates a data frame with time indices, labels, and steps for analysis functions
#'
#' @param start_idx Integer starting index for analysis
#' @param num_days Integer total number of days in dataset
#' @param time_labels Vector of time labels corresponding to each day
#' @return Data frame with columns: time_idx, time, t_step
create_time_grid <- function(start_idx, num_days, time_labels) {
  data.frame(
    time_idx = start_idx:num_days,
    time = time_labels[start_idx:num_days],
    t_step = start_idx:num_days
  )
}

#' Expand Time Grid for Multiple Pathogens
#'
#' Expands a time grid to include all pathogen-time combinations (sequential
#' pairing, not Cartesian product)
#'
#' @param time_grid Data frame from create_time_grid()
#' @param pathogen_names Character vector of pathogen names
#' @return Data frame with columns: pathogen, pathogen_idx, time_idx, time, t_step
expand_pathogen_grid <- function(time_grid, pathogen_names) {
  num_days <- nrow(time_grid)
  num_path <- length(pathogen_names)

  data.frame(
    pathogen = rep(pathogen_names, each = num_days),
    pathogen_idx = rep(seq_along(pathogen_names), each = num_days),
    time_idx = rep(time_grid$time_idx, times = num_path),
    time = rep(time_grid$time, times = num_path),
    t_step = rep(time_grid$t_step, times = num_path)
  )
}

#' Predict B-spline Basis Matrix
#'
#' Creates B-spline basis matrix for transforming spline coefficients to time series
#'
#' @param time Numeric vector of time points
#' @param days_per_knot Integer number of days between knots
#' @param spline_degree Integer degree of B-splines
#' @return Matrix B_true for transforming spline coefficients
#' @importFrom splines bs
#' @importFrom stats predict
predict_B_true <- function(time, days_per_knot, spline_degree) {
  X <- as.numeric(time)
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

#' Calculate Rt Denominator for Convolution
#'
#' Computes the denominator term for Rt calculations using generation interval convolution
#'
#' @param a Array of log-incidence posterior samples
#' @param i Integer time index for calculation
#' @param tau_max Integer maximum generation interval (days)
#' @param gi_dist Function returning generation interval probability for given day
#' @param is_multi_pathogen Logical indicating if this is a multi-pathogen model
#' @param pathogen_idx Integer pathogen index (NULL for total across pathogens)
#' @return Matrix of denominator values for Rt calculation
calc_rt_denominator <- function(a, i, tau_max, gi_dist, is_multi_pathogen = FALSE, pathogen_idx = NULL) {

  if (is_multi_pathogen && !is.null(pathogen_idx)) {
    # Multiple pathogen case - specific pathogen
    R_list <- matrix(0, nrow = dim(a)[1], ncol = 1)
    for(k in 0:(tau_max-1)) {
      R_list <- R_list + exp(a[, pathogen_idx, i-k]) * gi_dist(k)
    }
  } else if (is_multi_pathogen && is.null(pathogen_idx)) {
    # Multiple pathogen case - total across pathogens
    R_list <- matrix(0, nrow = dim(a)[1], ncol = 1)
    for(k in 0:(tau_max-1)) {
      R_list <- R_list + rowSums(exp(a[, , i-k])) * gi_dist(k)
    }
  } else {
    # Single pathogen case
    R_list <- matrix(0, nrow = length(a[,1]), ncol = 1)
    for(k in 0:(tau_max-1)) {
      R_list <- R_list + exp(a[, i-k]) * gi_dist(k)
    }
  }
  R_list
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
#' @param ... Additional arguments passed to calculation functions (e.g., tau_max, gi_dist for Rt)
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

  components <- extract_model_components(fitted_model)
  post <- rstan::extract(components$fit)

  # Transform data if using splines
  if (use_splines) {
    B_true <- predict_B_true(components$time, components$days_per_knot, components$spline_degree)
    a <- transform_posterior_multi(post, B_true, components$num_path, components$num_days)
  } else {
    a <- post$a
    # Ensure a is 3D array for multi-pathogen case
    if (length(dim(a)) == 2) {
      # If a is 2D, we need to reshape it to 3D
      # Assuming the structure is [samples, path*time] and needs to be [samples, path, time]
      a <- array(a, dim = c(nrow(a), components$num_path, components$num_days))
    }
  }

  time_grid <- create_time_grid(start_idx, components$num_days, components$time_labels)
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

  dplyr::bind_rows(pathogen_results, total_results) |>
    dplyr::arrange(pathogen != "Total", pathogen, t_step)
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

  components <- extract_model_components(fitted_model)
  post <- rstan::extract(components$fit)

  # Transform data if using splines
  if (use_splines) {
    B_true <- predict_B_true(components$time, components$days_per_knot, components$spline_degree)
    a <- transform_posterior_single(post, B_true, components$num_days)
  } else {
    a <- post$a
  }

  extra_args <- list(...)

  # Create results
  time_grid <- create_time_grid(start_idx, components$num_days, components$time_labels)

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

  return(results)
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
#' @param components Model components from extract_model_components()
#' @param extra_args List of additional arguments for calc_fn
#' @param threshold Numeric threshold for summary statistics
#' @return Data frame with expanded summary statistics
#' @importFrom purrr map2
calc_wrapper <- function (df, time_idx_col, pathogen_idx_col,
                          calc_fn, a, post, components, extra_args, threshold) {

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
