ps_get_a <- function (ps_fit,
                      X,
                      num_days = length(X),
                      time_labels,
                      days_per_knot = 5,
                      spline_degree = 3,
                      num_path = 4) {

  X <- as.numeric(X)

  knots <- get_knots(X, days_per_knot = days_per_knot,
                     spline_degree = spline_degree)
  num_knots <- length(knots)

  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)

  B_true <- splines::bs(
    seq(knots[1], knots[length(knots)], 1),
    knots = knots[2:(length(knots) - 1)],
    degree = spline_degree,
    intercept = TRUE
  )
  B_true <- t(predict(B_true, X))

  post <- rstan::extract(ps_fit)

  a <- array(data = NA, dim = c(nrow(post$a), num_path, num_days))

  for (j in 1:num_path) {
    print(i)
    for (k in 1:nrow(post$a)) {
      a[k, j, ] <- as.array(post$a[k, j, ]) %*% B_true
    }
  }

  a
}


ps_growth_rate_advantage <- function (a,
                                      X,
                                      num_days = length(X),
                                      time_labels,
                                      days_per_knot = 5,
                                      spline_degree = 3,
                                      num_path = 4,
                                      path_i = c(1, 2, 3),
                                      path_j = c(2, 3, 4),
                                      pathogen_names = c("Influenza A H3N2",
                                                         "Influenza A H1N1",
                                                         "Influenza B",
                                                         "Unknown")) {

  df <- data.frame()

  for (i in 2:num_days) {
    print(i)
    for (x in 1:length(path_i)) {
      j = path_i[x]
      k = path_j[x]

      quan <- quantile(
        (a[, j, i] - a[, j, (i - 1)]) -
          (a[, k, i] - a[, k, (i - 1)]),
        c(0.5, 0.025, 0.25, 0.75, 0.975)
      )
      prop <- length(
        a[, j, i][((a[, j, i] - a[, j, (i - 1)]) -
                     (a[, k, i] - a[, k, (i - 1)])) > 0]) /
        length(a[, j, i] - a[, j, (i - 1)])

      row <- data.frame(
        time = time_labels[i],
        t_step = i,
        y = quan[[1]],
        lb_50 = quan[[3]],
        ub_50 = quan[[4]],
        lb_95 = quan[[2]],
        ub_95 = quan[[5]],
        pathogen1 = pathogen_names[j],
        pathogen2 = pathogen_names[k],
        prop = prop
      )
      df <- rbind(df, row)

    }
  }
  df
}

ps_Rt_advantage <- function (a,
                             X,
                             num_days = length(X),
                             time_labels,
                             tau_max = 7,
                             gi_dist,
                             days_per_knot = 5,
                             spline_degree = 3,
                             num_path = 4,
                             path_i = c(1, 2, 3),
                             path_j = c(2, 3, 4),
                             pathogen_names = c("Influenza A H3N2",
                                                "Influenza A H1N1",
                                                "Influenza B",
                                                "Unknown")) {

  df <- data.frame()

  g_a <- sum(gi_dist(seq(0, tau_max - 1, 1)))

  for (i in tau_max:num_days) {
    print(i)

    for (x in 1:length(path_i)) {
      j1 <- path_i[x]
      j2 <- path_j[x]

      R_list1 <- matrix(0, nrow = length(a[, 1, 1]), ncol = 1)
      R_list2 <- matrix(0, nrow = length(a[, 1, 1]), ncol = 1)

      for (k in 0:(tau_max - 1)) {
        R_list1 <- R_list1 + exp(a[, j1, i - k]) * gi_dist(k)
        R_list2 <- R_list2 + exp(a[, j2, i - k]) * gi_dist(k)

      }
      R_list1 <- exp(a[, j1, i]) / (R_list1 / g_a)
      R_list2 <- exp(a[, j2, i]) / (R_list2 / g_a)

      quan <- quantile(R_list1 / R_list2, c(0.5, 0.025, 0.25, 0.75, 0.975))

      row <- data.frame(
        time = time_labels[i],
        t_step = i,
        y = quan[[1]],
        lb_50 = quan[[3]],
        ub_50 = quan[[4]],
        lb_95 = quan[[2]],
        ub_95 = quan[[5]],
        pathogen1 = pathogen_names[j1],
        pathogen2 = pathogen_names[j2]
      )
      df <- rbind(df, row)
    }
  }
  df
}

ps_Rt_advantage2 <- function (a,
                              X,
                              num_days = length(X),
                              time_labels,
                              tau_max = 7,
                              gi_dist_i,
                              gi_dist_j,
                              days_per_knot = 5,
                              spline_degree = 3,
                              num_path = 4,
                              path_i = c(1, 2, 3),
                              path_j = c(2, 3, 4),
                              pathogen_names = c("Influenza A H3N2",
                                                 "Influenza A H1N1",
                                                 "Influenza B",
                                                 "Unknown")) {

  df <- data.frame()

  for (i in tau_max:num_days) {
    print(i)

    for (x in 1:length(path_i)) {
      j1 <- path_i[x]
      j2 <- path_j[x]
      g_a1 <- sum(gi_dist_i[[x]](seq(0, tau_max - 1, 1)))
      g_a2 <- sum(gi_dist_j[[x]](seq(0, tau_max - 1, 1)))

      R_list1 <- matrix(0, nrow = length(a[, 1, 1]), ncol = 1)
      R_list2 <- matrix(0, nrow = length(a[, 1, 1]), ncol = 1)

      for (k in 0:(tau_max - 1)) {
        R_list1 <- R_list1 + exp(a[, j1, i - k]) * gi_dist_i[[x]](k)
        R_list2 <- R_list2 + exp(a[, j2, i - k]) * gi_dist_j[[x]](k)

      }
      R_list1 <- exp(a[, j1, i]) / (R_list1 / g_a1)
      R_list2 <- exp(a[, j2, i]) / (R_list2 / g_a2)

      quan <- quantile(R_list1 / R_list2, c(0.5, 0.025, 0.25, 0.75, 0.975))

      row <- data.frame(
        time = time_labels[i],
        t_step = i,
        y = quan[[1]],
        lb_50 = quan[[3]],
        ub_50 = quan[[4]],
        lb_95 = quan[[2]],
        ub_95 = quan[[5]],
        pathogen1 = pathogen_names[j1],
        pathogen2 = pathogen_names[j2]
      )
      df <- rbind(df, row)
    }
  }
  df
}
