#' Construct model
#'
#' @param method either random_walk() or p_spline()
#' @param pathogen_structure either single(), multiple(), or subtyped()
#' @param dow_effect logical whether to incorporate a day of week model
#'
#' @returns a list containing the data, the model parameters, and pathogen
#'  names
#' @export
#'
#' @examples
#' mod <- construct_model(
#'   method = random_walk(),
#'   pathogen_structure = single(
#'     case_timeseries = sarscov2$cases,
#'     time = sarscov2$date,
#'     pathogen_name = 'SARS-COV-2'
#'   ),
#'   dow_effect = TRUE
#' )
construct_model <- function (method,
                             pathogen_structure,
                             dow_effect = FALSE) {

  validate_class_inherits(
    method, 'EpiStrainDynamics.method'
  )
  validate_class_inherits(
    pathogen_structure, 'EpiStrainDynamics.pathogen_structure'
  )

  model_type <- get_model_type(
    method$method,
    pathogen_structure$pathogen_structure
  )

  time_seq <- seq(1, length(pathogen_structure$data$case_timeseries))

  # Set up day-of-week effects
  if (dow_effect) week_effect <- 7
  if (!dow_effect) week_effect <- 1
  DOW <- (time_seq %% week_effect) + 1

  model_params <- append(
    list(
      week_effect = week_effect,
      DOW = DOW
    ),
    pathogen_structure$model_params)

  if (method$method == 'p-spline') {
    knots <- get_knots(
      time_seq,
      days_per_knot = method$model_params$days_per_knot,
      spline_degree = method$model_params$spline_degree
    )
    model_params <- c(
      model_params,
      method$model_params,
      list(knots = knots)
    )
  }

  data <- append(
    list(time_seq = time_seq),
    pathogen_structure$data
  )

  model_input <- list(
    data = data,
    model_params = model_params,
    pathogen_names = pathogen_structure$pathogen_names,
    dow_effect = dow_effect
  )

  class(model_input) <- c(model_type, class(model_input))
  return(model_input)
}

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

#' Function for getting knot locations
#'
#' @param X days
#' @param days_per_knot days per knot
#' @param spline_degree spline degree
#'
#' @return location of knots
#'
get_knots <- function (X, days_per_knot = 3, spline_degree = 3) {

  X <- as.numeric(X)

  num_knots <- ceiling((max(X) - min(X)) / days_per_knot)

  first_knot <- min(X) - spline_degree * days_per_knot
  final_knot <- first_knot + days_per_knot * num_knots +
    2 * spline_degree * days_per_knot

  knots <- seq(first_knot, final_knot, by = days_per_knot)

  return(knots)
}
