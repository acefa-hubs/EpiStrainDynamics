
# Returns data.frame() of modeled growth rates
rw_single_growth_rate <- function(fit, num_days, time_labels){

  post <- rstan::extract(fit)

  df <- data.frame()

  for(i in 2:num_days){

    quan<- quantile(post$a[,i] - post$a[,(i-1)], c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(post$a[,i][(post$a[,i] - post$a[,(i-1)])>0])/length(post$a[,i] - post$a[,(i-1)])
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      prop = prop)
    df <- rbind(df, row)

  }

  df


}

# Returns data.frame() of modeled growth rates
rw_growth_rate <- function(fit, num_days, time_labels, num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){

  post <- rstan::extract(fit)

  df <- data.frame()

  for(i in 2:num_days){

    total_i   <- matrix(data=0, nrow=1, ncol=nrow(post$a))
    total_im1 <- matrix(data=0, nrow=1, ncol=nrow(post$a))

    for(j in 1:num_path){
      quan<- quantile(post$a[,j,i] - post$a[,j,(i-1)], c(0.5,0.025, 0.25, 0.75, 0.975))
      prop <- length(post$a[,j,i][(post$a[,j,i] - post$a[,j,(i-1)])>0])/length(post$a[,j,i] - post$a[,j,(i-1)])
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = pathogen_names[j],
                        prop = prop)
      df <- rbind(df, row)

      total_i   <- total_i   + exp(post$a[,j,i])
      total_im1 <- total_im1 + exp(post$a[,j,(i-1)])

    }

    quan<- quantile(log(total_i)-log(total_im1), c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(total_i[ (log(total_i)-log(total_im1)) >0])/length(log(total_i)-log(total_im1))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total",
                      prop = prop)
    df <- rbind(df, row)

  }

  df


}


# Returns data.frame() of modeled growth rates
ps_single_growth_rate <- function(fit,
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


  post <- rstan::extract(fit)

  a <- array(data=NA, dim=c(nrow(post$a), num_days))

  for(k in 1:nrow(post$a)){
    a[k,] <- as.array(post$a[k,]) %*% B_true
  }

  df <- data.frame()

  for(i in 2:num_days){
    quan<- quantile(a[,i] - a[,(i-1)], c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(a[,i][ (a[,i] - a[,(i-1)])>0])/length(a[,i] - a[,(i-1)])
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y=quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      prop = prop)
    df <- rbind(df, row)

  }

  df

}


# Returns data.frame() of modeled growth rates
ps_growth_rate <- function(fit,
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


  post <- rstan::extract(fit)

  a <- array(data=NA, dim=c(nrow(post$a), num_path, num_days))

  for(j in 1:num_path){
    for(k in 1:nrow(post$a)){
      a[k,j,] <- as.array(post$a[k,j,]) %*% B_true
    }
  }

  df <- data.frame()

  for(i in 2:num_days){

    total_i   <- matrix(data=0, nrow=1, ncol=nrow(a))
    total_im1 <- matrix(data=0, nrow=1, ncol=nrow(a))

    for(j in 1:num_path){
      quan<- quantile(a[,j,i] - a[,j,(i-1)], c(0.5,0.025, 0.25, 0.75, 0.975))
      prop <- length(a[,j,i][ (a[,j,i] - a[,j,(i-1)])>0])/length(a[,j,i] - a[,j,(i-1)])
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = pathogen_names[j],
                        prop = prop)
      df <- rbind(df, row)

      total_i   <- total_i   + exp(a[,j,i])
      total_im1 <- total_im1 + exp(a[,j,(i-1)])

    }

    quan<- quantile(log(total_i)-log(total_im1), c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(total_i[ (log(total_i)-log(total_im1)) >0])/length(log(total_i)-log(total_im1))
    row <- data.frame(time = time_labels[i],
                      t_step = i,
                      y = quan[[1]],
                      lb_50 = quan[[3]],
                      ub_50 = quan[[4]],
                      lb_95 = quan[[2]],
                      ub_95 = quan[[5]],
                      pathogen = "Total",
                      prop = prop)
    df <- rbind(df, row)

  }

  df



}


