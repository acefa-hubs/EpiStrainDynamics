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
