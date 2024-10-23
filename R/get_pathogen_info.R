get_pathogen_info <- function (data) {
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
  out <- list(pathogen_type = pathogen_type,
              pathogen_names = pathogen_names)
}
