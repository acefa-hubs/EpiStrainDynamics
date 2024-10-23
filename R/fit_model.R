fit_model <- function (data,
                       method = c('p-spline', 'random_walk'),
                       dow = TRUE,
                       covariance_structure = NULL,
                       iter,
                       warmup,
                       chains) {

  if ('influenzaA_subtypes' %in% names(data)) {
    pathogen_type <- 'influenzaA_subtypes'
    pathogen_names <- c(
      names(data$influenzaA_subtypes),
      names(data$main_pathogens)[names(data$main_pathogens) != 'influenzaA']
    )
  }
  if (!'influenzaA_subtypes' %in% names(data)) {
    pathogen_type <- 'other'
    pathogen_names <- names(data$main_pathogens)
  }

  # single pathogens
  if (length(pathogen_names) == 1 & method == 'p-spline') {
    model <- 'ps_single'
  }
  if (length(pathogen_names) == 1 & method == 'random_walk') {
    model <- 'rw_single'
  }

  # multiple pathogens
  if (length(pathogen_names) > 1) {

    if (method == 'p-spline' & pathogen_type == 'influenzaA_subtypes') {
      model <- 'ps_influenza'
    }
    if (method == 'p-spline' & pathogen_type == 'other') {
      model <- 'ps_mp'
    }

    if (method == 'random_walk' & pathogen_type == 'influenzaA_subtypes') {
      model <- 'rw_influenza'
    }
    if (method == 'random_walk' & pathogen_type == 'other') {
      model <- 'rw_mp'
    }

  }

  # data into format for stan
  Y <- data$cases
  main_pathogens <- t(matrix(
    data = c(
      data$main_pathogens$influenzaA,
      data$main_pathogens$influenzaB,
      data$main_pathogens$influenzaOther
    )
    ,
    ncol = 3
  ))
  influenzaA_subtypes <- t(matrix(
    data = c(
      data$influenzaA_subtypes$influenzaA.H3N2,
      data$influenzaA_subtypes$influenzaA.H1N1
    ),
    ncol = 2
  ))

  # fit stan model
  if (model == 'rw_influenza') {
    fit <- rw_influenza_stan(
      num_data = length(Y),
      num_path = length(pathogen_names),
      Y = Y,
      P1 = main_pathogens,
      P2 = influenzaA_subtypes,
      week_effect = 1,
      DOW = (df1$time %% 1) + 1,
      cov_structure = 1,
      iter = iter,
      warmup = warmup,
      chains = chains)
  }

  out <- list(fit = fit,
              model = model,
              pathogen_names = pathogen_names)

  class(out) <- 'EpiStrain.fit'
}
