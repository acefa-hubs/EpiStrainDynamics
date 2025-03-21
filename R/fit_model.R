fit_model <- function (constructed_model,
                       iter = 2000,
                       warmup = 1000,
                       chains = 3) {

  # if(!'EpiStrainDynamics.options' %in% class(options))
  #   stop("`options` must be created with model_options()")

  out <- fit(constructed_model, iter, warmup, chains)

  return(out)
}


fit <- function(constructed_model, iter, warmup, chains) UseMethod("fit")

# fit stan model
fit.rw_influenza <- function (constructed_model, iter, warmup, chains) {

  cases <- constructed_model$data$total_cases
  pathogen_names <- constructed_model$pathogen_names

  fit_object <- rw_influenza_stan(
    num_data = length(cases),
    num_path = length(pathogen_names),
    Y = cases,
    P1 = constructed_model$data$component_pathogens,
    P2 = constructed_model$data$influenzaA_subtypes,
    week_effect = constructed_model$model_params$week_effect,
    DOW = constructed_model$model_params$DOW,
    cov_structure = constructed_model$model_params$cov_structure,
    noise_structure = constructed_model$model_params$noise_structure,
    iter = iter,
    warmup = warmup,
    chains = chains)

  out <- list(fit = fit_object,
              constructed_model)

  class(out) <- 'EpiStrainDynamics.fit'

  return(out)
}
#
# fit.ps_influenza <- function () {
#   fit <- ps_influenza_stan(
#     num_data = length(cases),
#     num_path = length(pathogen_names),
#     num_knots = length(knots),
#     knots = knots,
#     spline_degree = spline_degree,
#     Y = cases,
#     X = time,
#     P1 = component_pathogens,
#     P2 = influenzaA_subtypes,
#     week_effect = week_effect,
#     DOW = DOW,
#     cov_structure = cov_structure,
#     noise_structure = noise_structure,
#     iter = iter,
#     warmup = warmup,
#     chains = chains)
# }
#
# if (model == 'rw_mp') {
#   fit <- rw_mp_stan(
#     num_data = length(cases),
#     num_path = length(pathogen_names),
#     Y = cases,
#     P = component_pathogens,
#     week_effect = week_effect,
#     DOW = DOW,
#     cov_structure = cov_structure,
#     noise_structure = noise_structure,
#     iter = iter,
#     warmup = warmup,
#     chains = chains
#   )
# }
# if (model == 'ps_mp') {
#   fit <- ps_mp_stan(
#     num_data = length(cases),
#     num_path = length(pathogen_names),
#     num_knots = length(knots),
#     knots = knots,
#     spline_degree = spline_degree,
#     Y = cases,
#     X = time,
#     P = component_pathogens,
#     week_effect = week_effect,
#     DOW = DOW,
#     cov_structure = cov_structure,
#     noise_structure = noise_structure,
#     iter = iter,
#     warmup = warmup,
#     chains = chains
#   )
# }
#
# if (model == 'ps_single') {
#   fit <- ps_single_stan(
#     num_data = length(cases),
#     num_knots = length(knots),
#     knots = knots,
#     spline_degree = spline_degree,
#     Y = cases,
#     X = time,
#     week_effect = week_effect,
#     DOW = DOW,
#     iter = iter,
#     warmup = warmup,
#     chains = chains
#   )
# }
# if (model == 'rw_single') {
#   fit <- rw_single_stan(
#     num_data = length(cases),
#     Y = cases,
#     week_effect = week_effect,
#     DOW = DOW
#   )
# }
