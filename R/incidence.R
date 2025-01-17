## Functions operating on random walk stan model fits

# Returns data.frame() of modeled incidence
rw_single_incidence <- function(rw_fit, num_days, time_labels){

  post <- rstan::extract(rw_fit)

  df <- data.frame()

  for(i in 1:num_days){
    quan <- quantile(post$a[,i], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=exp(quan[[1]]),
                      lb_50 = exp(quan[[3]]),
                      ub_50 = exp(quan[[4]]),
                      lb_95 = exp(quan[[2]]),
                      ub_95 = exp(quan[[5]]))
    df <- rbind(df, row)
  }

  df

}

# Returns data.frame() of modeled incidence
rw_incidence <- function(rw_fit, num_days, time_labels, num_path=4,
                         pathogen_names=c("Influenza A H3N2",
                                          "Influenza A H1N1",
                                          "Influenza B",
                                          "Unknown")){

  post <- rstan::extract(rw_fit)

  df <- data.frame()

  for(i in 1:num_days){

    total <- matrix(data=0, nrow=1, ncol=nrow(post$a))

    for(j in 1:num_path){
      quan<- quantile(post$a[,j,i], c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=exp(quan[[1]]),
                        lb_50 = exp(quan[[3]]),
                        ub_50 = exp(quan[[4]]),
                        lb_95 = exp(quan[[2]]),
                        ub_95 = exp(quan[[5]]),
                        pathogen = pathogen_names[j])
      df <- rbind(df, row)

      total <- total + exp(post$a[,j,i])

    }

    quan<- quantile(total, c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total")
    df <- rbind(df, row)

  }

  df


}


# Returns data.frame() of modeled incidence
ps_single_incidence <- function(ps_fit,
                                X,
                                num_days = length(X), time_labels,
                                days_per_knot = 5, spline_degree = 3){


  X <- as.numeric(X)

  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)


  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)

  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))


  post <- rstan::extract(ps_fit)

  a <- array(data=NA, dim=c(nrow(post$a), num_days))

  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }

  df <- data.frame()

  for(i in 1:num_days){
    total <- matrix(data=0, nrow=1, ncol=nrow(a))

    quan<- quantile(a[,i], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=exp(quan[[1]]),
                      lb_50 = exp(quan[[3]]),
                      ub_50 = exp(quan[[4]]),
                      lb_95 = exp(quan[[2]]),
                      ub_95 = exp(quan[[5]]))
    df <- rbind(df, row)

  }

  df


}


# Returns modelled incidence including day of the week effect
ps_single_incidence_dow <- function(ps_fit,
                                    X,
                                    num_days = length(X), time_labels,
                                    week_effect, DOW,
                                    days_per_knot = 5, spline_degree = 3){


  X <- as.numeric(X)

  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)


  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)

  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))


  post <- rstan::extract(ps_fit)

  a <- array(data=NA, dim=c(nrow(post$a), num_days))

  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }

  df <- data.frame()

  for(i in 1:num_days){
    total <- matrix(data=0, nrow=1, ncol=nrow(a))

    quan<- quantile(exp(a[,i])* week_effect*post$day_of_week_simplex[,DOW[i]], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]])
    df <- rbind(df, row)

  }

  df


}


# Returns data.frame() of modeled incidence
ps_incidence <- function(ps_fit,
                         X,
                         num_days = length(X), time_labels,
                         days_per_knot = 5, spline_degree = 3,
                         num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){


  X <- as.numeric(X)

  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)


  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)

  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))


  post <- rstan::extract(ps_fit)

  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))

  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }

  df <- data.frame()

  for(i in 1:num_days){

    total <- matrix(data=0, nrow=1, ncol=nrow(a))

    for(j in 1:num_path){
      quan<- quantile(a[,j,i], c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=exp(quan[[1]]),
                        lb_50 = exp(quan[[3]]),
                        ub_50 = exp(quan[[4]]),
                        lb_95 = exp(quan[[2]]),
                        ub_95 = exp(quan[[5]]),
                        pathogen = pathogen_names[j])
      df <- rbind(df, row)

      total <- total + exp(a[,j,i])

    }

    quan<- quantile(total, c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total")
    df <- rbind(df, row)

  }

  df


}



ps_incidence_dow <- function(ps_fit,
                             X,
                             num_days = length(X), time_labels, DOW, week_effect,
                             days_per_knot = 5, spline_degree = 3,
                             num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){


  X <- as.numeric(X)

  knots <- get_knots(X, days_per_knot = days_per_knot, spline_degree = spline_degree)
  num_knots <- length(knots)


  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)

  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1), knots = knots[2:(length(knots)-1)], degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))


  post <- rstan::extract(ps_fit)

  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))

  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }

  df <- data.frame()

  for(i in 1:num_days){

    total <- matrix(data=0, nrow=1, ncol=nrow(a))

    for(j in 1:num_path){
      quan<- quantile(exp(a[,j,i])* week_effect*post$day_of_week_simplex[,DOW[i]], c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = pathogen_names[j])
      df <- rbind(df, row)

      total <- total + exp(a[,j,i])

    }

    quan<- quantile(total*week_effect*post$day_of_week_simplex[,DOW[i]], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total")
    df <- rbind(df, row)

  }

  df


}






