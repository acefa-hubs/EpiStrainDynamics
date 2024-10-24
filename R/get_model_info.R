get_model_info <- function (pathogen_type, pathogen_names, method) {

  # single pathogens
  if (is.null(pathogen_names) & method == 'p-spline') {
    model <- 'ps_single'
  }
  if (is.null(pathogen_names) & method == 'random_walk') {
    model <- 'rw_single'
  }

  # multiple pathogens
  if (!is.null(pathogen_names)) {

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
  return(model)
}
