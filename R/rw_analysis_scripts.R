## Functions operating on random walk stan model fits

# Returns data.frame() of modeled incidence
rw_incidence <- function(rw_fit, num_days, time_labels, num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){
  
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


# Returns data.frame() of modeled relative proportions of different pathogens (for comparison to data)
rw_proportion <- function(rw_fit,
                          num_days, time_labels,
                          num_path=4,
                          comb_num=list(c(1,2,3), c(1,2), c(3), c(1), c(2)),
                          comb_den=list(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4), c(1,2), c(1,2) ),
                          comb_names=c("Influenza", "Influenza A", "Influenza B", "H3N2", "H1N1")){
  
  post <- rstan::extract(rw_fit)
  
  df <- data.frame()
  
  for(i in 1:num_days){
    
    total <- matrix(data=0, nrow=1, ncol=nrow(post$a))
    
    for(j in 1:length(comb_num)){
      num_index <- comb_num[[j]]
      den_index <- comb_den[[j]]
      
      if(length(num_index)>1){
        num <- rowSums(exp(post$a[,num_index,i]))
      } else{
        num <- exp(post$a[,num_index,i])
      }
      
      den <- rowSums(exp(post$a[,den_index,i]))
      
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
rw_growth_rate <- function(rw_fit, num_days, time_labels, num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){
  
  post <- rstan::extract(rw_fit)
  
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


# Modeled Rt estimates from incidence curve
rw_Rt <- function(rw_fit,
                  num_days, time_labels,
                  tau_max=7,
                  gi_dist,
                  num_path=4, pathogen_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B", "Unknown")){
  
  post <- rstan::extract(rw_fit)
  
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
