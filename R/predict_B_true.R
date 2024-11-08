predict_B_true <- function (time, days_per_knot, spline_degree) {

  time <- as.numeric(time)

  knots <- get_knots(
    time,
    days_per_knot = days_per_knot,
    spline_degree = spline_degree
  )
  num_knots <- length(knots)
  num_basis <- num_knots + spline_degree - 1

  B_true <- splines::bs(
    seq(knots[1], knots[length(knots)], 1),
    knots = knots[2:(length(knots) - 1)],
    degree = spline_degree,
    intercept = TRUE
  )
  B_true <- t(predict(B_true, time))

  return(B_true)
}
