
## Functions operating on random walk stan model fits

# Function for getting knot locations
get_knots <- function(X, days_per_knot, spline_degree=3){
   
  X <- as.numeric(X)
  
  num_knots <- ceiling((max(X)-min(X))/days_per_knot)
  
  first_knot <- min(X) - spline_degree*days_per_knot
  final_knot <- first_knot + days_per_knot*num_knots + 2*spline_degree*days_per_knot
  
  knots <- seq(first_knot, final_knot, by=days_per_knot)
  
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






# Returns data.frame() of modeled relative proportions of different pathogens (for comparison to data)
ps_proportion <- function(ps_fit,
                         X,
                         num_days = length(X), time_labels,
                         days_per_knot = 5, spline_degree = 3,
                         num_path=4, 
                         comb_num=list(c(1,2,3), c(1,2), c(3), c(1), c(2)),
                         comb_den=list(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4), c(1,2), c(1,2) ),
                         comb_names=c("Influenza", "Influenza A", "Influenza B", "H3N2", "H1N1") ){
  
  
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
    
    for(j in 1:length(comb_num)){
      num_index <- comb_num[[j]]
      den_index <- comb_den[[j]]
      
      if(length(num_index)>1){
        num <- rowSums(exp(a[,num_index,i]))
      } else{
        num <- exp(a[,num_index,i])
      }
      
      den <- rowSums(exp(a[,den_index,i]))
      
      quan<- quantile(num/den, c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = time_labels[i],
                        t_step = i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = comb_names[j])
      df <- rbind(df, row)
      
    }
    
  }
  
  df
  
  
}



# Returns data.frame() of modeled growth rates
ps_growth_rate <- function(ps_fit,
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
ps_Rt <- function(ps_fit,
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
  
  
  post <- rstan::extract(ps_fit)
  
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



#####################################################################################
#####################################################################################

ps_growth_rate_advantage <- function(ps_fit,
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
  
  for(i in 2:num_days){
    print(i)
    for(j in 1:num_path){
      for(k in 1:num_path){
        
        quan<- quantile( (a[,j,i] - a[,j,(i-1)])-(a[,k,i] - a[,k,(i-1)]), c(0.5,0.025, 0.25, 0.75, 0.975))
        prop <- length(a[,j,i][ ( (a[,j,i] - a[,j,(i-1)]) - (a[,k,i] - a[,k,(i-1)]) )>0])/length(a[,j,i] - a[,j,(i-1)])
        row <- data.frame(time = time_labels[i],
                          t_step = i,
                          y=quan[[1]],
                          lb_50 = quan[[3]],
                          ub_50 = quan[[4]],
                          lb_95 = quan[[2]],
                          ub_95 = quan[[5]],
                          pathogen1 = pathogen_names[j],
                          pathogen2 = pathogen_names[k],
                          prop = prop)
        df <- rbind(df, row)
        
      }
      
    }
    
  }
  
  df
  
}







ps_forecast <- function(ps_fit,
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
  
  total <- matrix(data=0, ncol=28, nrow=nrow(a))
  for(j in 1:num_path){
    gr <- (a[,j,num_days] - a[,j,(num_days-7)])/7
    for_inc <- exp(a[,j,num_days] + gr%*%t(seq(1,28,1)))
    total<- total+for_inc
    
    for(i in 1:28){
      
      quan<- quantile(for_inc[,i], c(0.5,0.025, 0.25, 0.75, 0.975))
      row <- data.frame(time = i +time_labels[num_days],
                        step= i,
                        y=quan[[1]],
                        lb_50 = quan[[3]],
                        ub_50 = quan[[4]],
                        lb_95 = quan[[2]],
                        ub_95 = quan[[5]],
                        pathogen = pathogen_names[j])
      df <- rbind(df, row)
      
    }
    
    
  }
  
  for(i in 1:28){
    
    quan<- quantile(total[,i], c(0.5,0.025, 0.25, 0.75, 0.975))
    row <- data.frame(time = i +time_labels[num_days],
                      step=i,
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












