## Functions operating on random walk stan model fits


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

