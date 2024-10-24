fit_model <- function (data,
                       method = c('p-spline', 'random_walk'),
                       spline_degree = NULL,
                       days_per_knot = NULL,
                       covariance_structure = NULL,
                       iter,
                       warmup,
                       chains) {

  try(if (model == 'p-spline' & is.null(spline_degree))
    stop("Must specify spline_degree if using p-spline method"))

  pathogen_type <- get_pathogen_info(data)$pathogen_type
  pathogen_names <- get_pathogen_info(data)$pathogen_names
  model <- get_model_info(pathogen_type, pathogen_names, method)

  # data into format for stan
  cases <- data$cases
  time <- data$time
  if ('main_pathogens' %in% names(data)) {
    main_pathogens <- t(matrix(
      data = do.call(c, data$main_pathogens),
      ncol = length(data$main_pathogens)
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

  # fit stan model
  if (model == 'rw_influenza') {
    fit <- rw_influenza_stan(
      num_data = length(cases),
      num_path = length(pathogen_names),
      Y = cases,
      P1 = main_pathogens,
      P2 = influenzaA_subtypes,
      week_effect = 1,
      DOW = (time %% 1) + 1,
      cov_structure = 1,
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
      P1 = main_pathogens,
      P2 = influenzaA_subtypes,
      week_effect = 1,
      DOW = (time %% 1) + 1,
      cov_structure = 1,
      iter = iter,
      warmup = warmup,
      chains = chains)
  }

  if (model == 'rw_mp') {
    fit <- rw_mp_stan(
      num_data = length(cases),
      num_path = length(pathogen_names),
      Y = cases,
      P = main_pathogens,
      week_effect = 1,
      DOW = (time %% 1) + 1,
      cov_structure = 1,
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
      P = main_pathogens,
      week_effect = 1,
      DOW = (time %% 1) + 1,
      cov_structure = 0,
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
      week_effect = 7,
      DOW = (time %% 7) + 1,
      iter = iter,
      warmup = warmup,
      chains = chains
    )
  }
  if (model == 'rw_single') {
    fit <- rw_single_stan(
      num_data = length(cases),
      Y = cases,
      week_effect = 1,
      DOW = (time %% 1) + 1
    )
  }

  out <- list(fit = fit,
              model = model,
              pathogen_names = pathogen_names)
  if (exists(knots)) {
    out <- c(out, knots = knots)
  }

  class(out) <- 'EpiStrain.fit'
  return(out)
}
