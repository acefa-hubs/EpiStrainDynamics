# Rt functions
# Modeled Rt estimates from incidence curve
rw_Rt <- function(fit,
                  num_days, time_labels,
                  tau_max=7,
                  gi_dist,
                  num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){

  post <- rstan::extract(fit)

  df <- data.frame()

  g_a <- sum(gi_dist(seq(0,tau_max-1,1)))

  for(i in tau_max:num_days){

    R_list_T <- matrix(0, nrow=length(post$a[,1,1]),ncol=1) #
    for(j in 1:num_path){

      R_list <- matrix(0, nrow=length(post$a[,1,1]),ncol=1)

      for(k in 0:(tau_max-1)){
        R_list <- R_list + exp(post$a[,j,i-k])*gi_dist(k)

        R_list_T <- R_list_T + exp(post$a[,j,i-k])*gi_dist(k)#
      }
      R_list <- exp(post$a[,j,i])/(R_list/g_a)

      quan<- quantile(R_list, c(0.5,0.025, 0.25, 0.75, 0.975))
      prop <- length(R_list[R_list>1]) / length(R_list)
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

    }


    R_list_T <- rowSums(exp(post$a[ , ,i])) /(R_list_T/g_a)
    quan<- quantile(R_list_T, c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(R_list_T[R_list_T>1]) / length(R_list_T)
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



# Modeled Rt estimates from incidence curve
rw_single_Rt <- function(fit,
                         num_days, time_labels,
                         tau_max=7,
                         gi_dist){

  post <- rstan::extract(fit)

  df <- data.frame()

  g_a <- sum(gi_dist(seq(0,tau_max-1,1)))

  for(i in tau_max:num_days){

    R_list <- matrix(0, nrow=length(post$a[,1]),ncol=1)

    for(k in 0:(tau_max-1)){
      R_list <- R_list + exp(post$a[,i-k])*gi_dist(k)
    }
    R_list <- exp(post$a[,i])/(R_list/g_a)

    quan<- quantile(R_list, c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(R_list[R_list>1]) / length(R_list)
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
ps_single_Rt <- function(fit,
                         X,
                         num_days = length(X), time_labels,
                         tau_max = 7,
                         gi_dist,
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

  g_a <- sum(gi_dist(seq(0,tau_max-1,1)))

  for(i in tau_max:num_days){

    R_list <- matrix(0, nrow=length(a[,1]),ncol=1)

    for(k in 0:(tau_max-1)){
      R_list <- R_list + exp(a[,i-k])*gi_dist(k)
    }
    R_list <- exp(a[,i])/(R_list/g_a)

    quan<- quantile(R_list, c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(R_list[R_list>1]) / length(R_list)
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
ps_Rt <- function(fit,
                  X,
                  num_days = length(X), time_labels,
                  tau_max = 7,
                  gi_dist,
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

  g_a <- sum(gi_dist(seq(0,tau_max-1,1)))

  for(i in tau_max:num_days){

    R_list_T <- matrix(0, nrow=length(a[,1,1]),ncol=1) #
    for(j in 1:num_path){

      R_list <- matrix(0, nrow=length(a[,1,1]),ncol=1)

      for(k in 0:(tau_max-1)){
        R_list <- R_list + exp(a[,j,i-k])*gi_dist(k)

        R_list_T <- R_list_T + exp(a[,j,i-k])*gi_dist(k)#
      }
      R_list <- exp(a[,j,i])/(R_list/g_a)

      quan<- quantile(R_list, c(0.5,0.025, 0.25, 0.75, 0.975))
      prop <- length(R_list[R_list>1]) / length(R_list)
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

    }

    R_list_T <- rowSums(exp(a[ , ,i])) /(R_list_T/g_a)
    quan<- quantile(R_list_T, c(0.5,0.025, 0.25, 0.75, 0.975))
    prop <- length(R_list_T[R_list_T>1]) / length(R_list_T)
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

