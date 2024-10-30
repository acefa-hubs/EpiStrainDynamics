# Defining functions needed for producing SFig 14
###############################################################################

ps_get_a <- function(post,
                     X,
                     num_days = length(X),
                     time_labels,
                     min_time,
                     max_time,
                     days_per_knot = 5,
                     spline_degree = 3,
                     num_path = 4,
                     pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")) {
  X <- as.numeric(X)

  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
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

  a <- array(data = NA, dim = c(nrow(post$a), num_path, num_days))

  for (j in 1:num_path) {
    print(j)
    for (k in 1:nrow(post$a)) {
      a[k, j, ] <- as.array(post$a[k, j, ]) %*% B_true
    }
  }

  a

}


ps_cumulative_incidence <- function(a,
                                    X,
                                    num_days = length(X),
                                    time_labels,
                                    min_time,
                                    max_time,
                                    days_per_knot = 5,
                                    spline_degree = 3,
                                    num_path = 4,
                                    pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")) {
  df <- data.frame()

  min_num_day <- which(time_labels == min_time)
  max_num_day <- which(time_labels == max_time)

  total <- rowSums(exp(a[, 1, min_num_day:max_num_day])) +
    rowSums(exp(a[, 2, min_num_day:max_num_day])) +
    rowSums(exp(a[, 3, min_num_day:max_num_day])) +
    rowSums(exp(a[, 4, min_num_day:max_num_day]))
  for (j in 1:num_path) {
    quan <- quantile(rowSums(exp(a[, j, min_num_day:max_num_day])), c(0.5, 0.025, 0.25, 0.75, 0.975))

    row <- data.frame(
      time_min = min_time,
      time_max = max_time,
      y = quan[[1]],
      lb_50 = quan[[3]],
      ub_50 = quan[[4]],
      lb_95 = quan[[2]],
      ub_95 = quan[[5]],
      pathogen = pathogen_names[j],
      label = "total"
    )
    df <- rbind(df, row)

    quan <- quantile((rowSums(exp(a[, j, min_num_day:max_num_day])) / total), c(0.5, 0.025, 0.25, 0.75, 0.975))

    row <- data.frame(
      time_min = min_time,
      time_max = max_time,
      y = quan[[1]],
      lb_50 = quan[[3]],
      ub_50 = quan[[4]],
      lb_95 = quan[[2]],
      ub_95 = quan[[5]],
      pathogen = pathogen_names[j],
      label = "proportion"
    )
    df <- rbind(df, row)

  }

  df

}

raw_cumulative_incidence <- function(df, min_time, max_time) {
  df_tmp <- df[df$date >= min_time & df$date <= max_time, ]

  sero1 <- sum(df_tmp$sero1)
  sero2 <- sum(df_tmp$sero2)
  sero3 <- sum(df_tmp$sero3)
  sero4 <- sum(df_tmp$sero4)
  total <- sero1 + sero2 + sero3 + sero4

  mp <- MultinomCI(c(sero1, sero2, sero3, sero4))

  # propCI
  df_r <- data.frame()
  for (i in 1:4) {
    row_df <- data.frame(
      time_min = min_time,
      time_max = max_time,
      y = mp[i, 1],
      lb_95 = mp[i, 2],
      ub_95 = mp[i, 3],
      lb_50 = mp[i, 2],
      ub_50 = mp[i, 3],
      pathogen = as.character(i),
      label = "proportion"
    )

    df_r <- rbind(df_r, row_df)
  }
  df_r

}
