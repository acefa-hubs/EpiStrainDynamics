fit_model <- function (data,
                       method = c(
                         'p-spline',
                         'random_walk'),
                       spline_degree = NULL,
                       days_per_knot = NULL,
                       dow_effect = TRUE,
                       smoothing_structure = c(
                         'single',
                         'independent',
                         'correlated'),
                       iter,
                       warmup,
                       chains) {

  # add this check for days_per_knot as well
  pspline_dependencies <- method == 'p-spline' &
    (is.null(spline_degree) | is.null(days_per_knot))
  try (if (pspline_dependencies)
    stop("Must specify spline_degree if using p-spline method"))

  pathogen_type <- get_pathogen_info(data)$pathogen_type
  pathogen_names <- get_pathogen_info(data)$pathogen_names
  model <- get_model_info(pathogen_type, pathogen_names, method)

  # data into format for stan
  cases <- data$cases
  # need to set checks or internal code to modify format so that it's integer
  time <- data$time
  if ('component_pathogens' %in% names(data)) {
    component_pathogens <- t(matrix(
      data = do.call(c, data$component_pathogens),
      ncol = length(data$component_pathogens)
    ))
  }
  if ('influenzaA_subtypes' %in% names(data)) {
    influenzaA_subtypes <- t(matrix(
      data = do.call(c, data$influenzaA_subtypes),
      ncol = length(data$influenzaA_subtypes)
    ))
  }

  if (method == 'p-spline') {
    # Calculate the locations of equally spaced knots
    knots <- get_knots(time, days_per_knot = days_per_knot,
                       spline_degree = spline_degree)
  }

  if (dow_effect) {
    week_effect <- 7
    DOW <- (time %% 7) + 1
  }
  if (!dow_effect) {
    week_effect <- 1
    DOW <- (time %% 1) + 1
  }

  if (smoothing_structure == 'single') cov_structure <- 0
  if (smoothing_structure == 'independent') cov_structure <- 1
  if (smoothing_structure == 'correlated') cov_structure <- 2

  # fit stan model
  if (model == 'rw_influenza') {
    fit <- rw_influenza_stan(
      num_data = length(cases),
      num_path = length(pathogen_names),
      Y = cases,
      P1 = component_pathogens,
      P2 = influenzaA_subtypes,
      week_effect = week_effect,
      DOW = DOW,
      cov_structure = cov_structure,
      iter = iter,
      warmup = warmup,
      chains = chains)
  }
  if (model == 'ps_influenza') {
    fit <- ps_influenza_stan(
      num_data = length(cases),
      num_path = length(pathogen_names),
      num_knots = length(knots),
      knots = knots,
      spline_degree = spline_degree,
      Y = cases,
      X = time,
      P1 = component_pathogens,
      P2 = influenzaA_subtypes,
      week_effect = week_effect,
      DOW = DOW,
      cov_structure = cov_structure,
      iter = iter,
      warmup = warmup,
      chains = chains)
  }

  if (model == 'rw_mp') {
    fit <- rw_mp_stan(
      num_data = length(cases),
      num_path = length(pathogen_names),
      Y = cases,
      P = component_pathogens,
      week_effect = week_effect,
      DOW = DOW,
      cov_structure = cov_structure,
      iter = iter,
      warmup = warmup,
      chains = chains
    )
  }
  if (model == 'ps_mp') {
    fit <- ps_mp_stan(
      num_data = length(cases),
      num_path = length(pathogen_names),
      num_knots = length(knots),
      knots = knots,
      spline_degree = spline_degree,
      Y = cases,
      X = time,
      P = component_pathogens,
      week_effect = week_effect,
      DOW = DOW,
      cov_structure = cov_structure,
      iter = iter,
      warmup = warmup,
      chains = chains
    )
  }

  if (model == 'ps_single') {
    fit <- ps_single_stan(
      num_data = length(cases),
      num_knots = length(knots),
      knots = knots,
      spline_degree = spline_degree,
      Y = cases,
      X = time,
      week_effect = week_effect,
      DOW = DOW,
      iter = iter,
      warmup = warmup,
      chains = chains
    )
  }
  if (model == 'rw_single') {
    fit <- rw_single_stan(
      num_data = length(cases),
      Y = cases,
      week_effect = week_effect,
      DOW = DOW
    )
  }

  out <- list(fit = fit,
              model = model,
              pathogen_names = pathogen_names)
  if (exists('knots')) {
    out <- c(out, knots = knots)
  }

  class(out) <- 'EpiStrain.fit'
  return(out)
}
