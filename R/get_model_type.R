#' Extract model type from method and pathogen structure
#'
#' @param method method function, either random_walk() or p_spline()
#' @param pathogens pathogen structure, either single, multiple, or subtyped
#'
#' @returns model type
#'
get_model_type <- function (method, pathogens) {

  is_ps_single <- method == 'p-spline' & pathogens == 'single'
  is_rw_single <- method == 'random-walk' & pathogens == 'single'
  is_ps_multiple <- method == 'p-spline' & pathogens == 'multiple'
  is_rw_multiple <- method == 'random-walk' & pathogens == 'multiple'
  is_ps_subtyped <- method == 'p-spline' & pathogens == 'subtyped'
  is_rw_subtyped <- method == 'random-walk' & pathogens == 'subtyped'

  if (is_ps_single) model_type <- 'ps_single'
  if (is_rw_single) model_type <- 'rw_single'
  if (is_ps_multiple) model_type <- 'ps_multiple'
  if (is_rw_multiple) model_type <- 'rw_multiple'
  if (is_ps_subtyped) model_type <- 'ps_subtyped'
  if (is_rw_subtyped) model_type <- 'rw_subtyped'

  return(model_type)
}

#' Get covariate structure from assigned smoothing structure
#'
#' @param smoothing_structure either 'shared' (all pathogens have the same;
#'   tau[1]), 'independent' (each pathogen has completely independent smoothing
#'   structure; tau[number of pathogens]), or 'correlated' (smoothing structure is
#'   correlated among pathogens Sigma[number of pathogens, number of pathogens])
#'
#' @returns numeric value for covariance structure needed by stan models
#'
get_cov_structure <- function (smoothing_structure) {
  if (smoothing_structure == 'single') cov_structure <- 0
  if (smoothing_structure == 'independent') cov_structure <- 1
  if (smoothing_structure == 'correlated') cov_structure <- 2
  return(cov_structure)
}

#' Get noise structure from specified observation noise
#'
#' @param observation_noise either 'observation_noise_only' (only includes
#'   observation noise - the same between pathogens) or pathogen_specific_noise'
#'   (includes noise in individual pathogens as well)
#'
#' @returns numeric value for noise structure needed by stan models
#'
get_noise_structure <- function (observation_noise) {
  if (observation_noise == 'observation_noise_only') noise_structure <- 0
  if (observation_noise == 'pathogen_specific_noise') noise_structure <- 1
  return(noise_structure)
}
