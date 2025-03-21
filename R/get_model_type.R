get_model_type <- function (method, pathogens) {

  is_ps_single <- method == 'p-spline' & pathogens == 'single'
  is_rw_single <- method == 'random-walk' & pathogens == 'single'
  is_ps_mp <- method == 'p-spline' & pathogens == 'multiple'
  is_rw_mp <- method == 'random-walk' & pathogens == 'multiple'
  is_ps_influenza <- method == 'p-spline' & pathogens == 'subtyped'
  is_rw_influenza <- method == 'random-walk' & pathogens == 'subtyped'

  if (is_ps_single) model_type <- 'ps_single'
  if (is_rw_single) model_type <- 'rw_single'
  if (is_ps_mp) model_type <- 'ps_mp'
  if (is_rw_mp) model_type <- 'rw_mp'
  if (is_ps_influenza) model_type <- 'ps_influenza'
  if (is_rw_influenza) model_type <- 'rw_influenza'

  return(model_type)
}

get_cov_structure <- function (smoothing_structure) {
  if (smoothing_structure == 'single') cov_structure <- 0
  if (smoothing_structure == 'independent') cov_structure <- 1
  if (smoothing_structure == 'correlated') cov_structure <- 2
  return(cov_structure)
}

get_noise_structure <- function (observation_noise) {
  if (observation_noise == 'observation_noise_only') noise_structure <- 0
  if (observation_noise == 'pathogen_specific_noise') noise_structure <- 1
  return(noise_structure)
}
