get_model_info <- function (pathogen_type, pathogen_names, method) {

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
  return(model)
}
