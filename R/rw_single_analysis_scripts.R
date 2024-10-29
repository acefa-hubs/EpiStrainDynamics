## Functions operating on random walk stan model fits

# Returns data.frame() of modeled incidence
rw_single_incidence <- function(rw_fit, num_days, time_labels){
  
  post <- rstan::extract(rw_fit)
  
  df <- data.frame()
  
  for(i in 1:num_days){
    quan<- quantile(post$a[,i], c(0.5,0.025, 0.25, 0.75, 0.975))
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


# Returns data.frame() of modeled growth rates
rw_single_growth_rate <- function(rw_fit, num_days, time_labels){
  
  post <- rstan::extract(rw_fit)
  
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


# Modeled Rt estimates from incidence curve
rw_Rt <- function(rw_fit,
                  num_days, time_labels,
                  tau_max=7,
                  gi_dist){
  
  post <- rstan::extract(rw_fit)
  
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
