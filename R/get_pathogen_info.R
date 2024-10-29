get_pathogen_info <- function (data) {
  if ('influenzaA_subtypes' %in% names(data)) {
    pathogen_type <- 'influenzaA_subtypes'
    message('Ensure that component influenza A is listed first.')
    pathogen_names <- c(
      names(data$influenzaA_subtypes),
      names(data$component_pathogens)[-1]
    )
  }
  if (!'influenzaA_subtypes' %in% names(data)) {
    pathogen_type <- 'other'
    pathogen_names <- names(data$component_pathogens)
  }
  out <- list(pathogen_type = pathogen_type,
              pathogen_names = pathogen_names)
}
