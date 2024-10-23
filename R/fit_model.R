fit_model <- function (data,
                       method = c('p-spline', 'random_walk'),
                       dow = TRUE,
                       covariance_structure = NULL,
                       iter,
                       warmup,
                       chains) {

  pathogen_type <- get_pathogen_info(data)$pathogen_type
  pathogen_names <- get_pathogen_info(data)$pathogen_names
  model <- get_model_info(pathogen_type, pathogen_names, method)

  # data into format for stan
  cases <- data$cases
  time <- data$time
  main_pathogens <- t(matrix(data = do.call(c, data$main_pathogens), ncol = 3))
  if (pathogen_type == 'influenzaA_subtypes') {
    influenzaA_subtypes <- t(matrix(
      data = do.call(c, data$influenzaA_subtypes),
      ncol = 2
    ))
  }

  if (method == 'p-spline') {
    # Calculate the locations of equally spaced knots
    knots <- get_knots(time, days_per_knot = 5, spline_degree = 3)
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
  if (model == 'ps_mp') {
    fit <- ps_mp_stan(
      num_data = length(cases),
      num_path = length(pathogen_names),
      num_knots = length(knots),
      knots = knots,
      spline_degree = 3,
      Y = cases,
      X = time,
      P = main_pathogens,
      week_effect = 1,
      DOW = (time %% 1) + 1,
      cov_structure = 0,
      iter = iter,
      warmup = warmup,
      chains = chains)
  }
  if (model == 'ps_single') {
    ps_single_fit <- ps_single_stan(
      num_data = length(cases),
      num_knots = length(knots),
      knots = knots,
      spline_degree = 3,
      Y = cases,
      X = time,
      week_effect = 7,
      DOW = (time %% 7) + 1,
      iter = iter,
      warmup = warmup,
      chains = chains)
  }

  out <- list(fit = fit,
              model = model,
              pathogen_names = pathogen_names)

  class(out) <- 'EpiStrain.fit'
}
