### Loading required packages
library("MASS")
library("randomForest")
library("data.table")
library("permimp")

source("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/empirical_weight.R")
library(doParallel)
library(foreach)
registerDoParallel(80)
t_95 = qt(.975, 2000-1)

rho = .9

sim_100_withforest = foreach(s = 1:160, .combine = "rbind")%dopar%{
  # TT1 = Sys.time()
  set.seed(s+2023)
  
  n    = 2000
  rho  = 0.9
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X = cbind(X_12, X_310)
  colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y = c(X%*%beta + rnorm(n, 0, 1))
  
  reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  full_model    = randomForest(X[,], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  
  set.seed(s+347)
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X_2 = cbind(X_12, X_310)
  colnames(X_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y_2 = c(X_2%*%beta + rnorm(n, 0, 1))
  
  
  N = 200
  zero_term = c()
  first_term = c()
  second_term = c()
  third_term = c()
  for (i in 1:n){
    beta_curr = beta[-1]
    beta_curr[1] = (1+rho)*beta_curr[1]
    curr_y = rnorm(N, beta_curr%*%X_2[i,-1],sqrt(2-rho^2))
    curr_x = rnorm(N, X_2[i,2]*rho, sqrt(1-rho^2))
    
    second_data =  as.data.frame(cbind(curr_x, matrix(rep(X_2[i,-1], N), ncol = 9, byrow = T)))
    colnames(second_data) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
    
    zero_term[i] = 2*mean((curr_y - predict(full_model, X_2[i,]))*
                            (Y_2[i] - predict(full_model, X_2[i,])))
    first_term[i] = mean((curr_y  - predict(full_model, X_2[i,]))^2)
    second_term[i] = mean((curr_y -  predict(full_model, newdata = second_data))^2)
    third_term[i] = mean((Y_2[i] -  predict(full_model, newdata = second_data))^2)
  }
  #### Variance and Std
  base_perm = mean( first_term - second_term + third_term)
  var_est = mean(( first_term - second_term + third_term-base_perm)^2)
  c(mean(zero_term), mean(first_term), mean(second_term), mean(third_term), var_est)
}

write.csv(sim_100_withforest, "~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-15/var_explore_2_1_nof_rf_large_rho9.csv" )



rho = .7

sim_100_withforest = foreach(s = 1:160, .combine = "rbind")%dopar%{
  # TT1 = Sys.time()
  set.seed(s+2023)
  
  n    = 2000
  rho  = 0.7
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X = cbind(X_12, X_310)
  colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y = c(X%*%beta + rnorm(n, 0, 1))
  
  reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  full_model    = randomForest(X[,], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  
  set.seed(s+347)
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X_2 = cbind(X_12, X_310)
  colnames(X_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y_2 = c(X_2%*%beta + rnorm(n, 0, 1))
  
  
  N = 200
  zero_term = c()
  first_term = c()
  second_term = c()
  third_term = c()
  for (i in 1:n){
    beta_curr = beta[-1]
    beta_curr[1] = (1+rho)*beta_curr[1]
    curr_y = rnorm(N, beta_curr%*%X_2[i,-1],sqrt(2-rho^2))
    curr_x = rnorm(N, X_2[i,2]*rho, sqrt(1-rho^2))
    
    second_data =  as.data.frame(cbind(curr_x, matrix(rep(X_2[i,-1], N), ncol = 9, byrow = T)))
    colnames(second_data) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
    
    zero_term[i] = 2*mean((curr_y - predict(full_model, X_2[i,]))*
                            (Y_2[i] - predict(full_model, X_2[i,])))
    first_term[i] = mean((curr_y  - predict(full_model, X_2[i,]))^2)
    second_term[i] = mean((curr_y -  predict(full_model, newdata = second_data))^2)
    third_term[i] = mean((Y_2[i] -  predict(full_model, newdata = second_data))^2)
  }
  #### Variance and Std
  base_perm = mean( first_term - second_term + third_term)
  var_est = mean(( first_term - second_term + third_term-base_perm)^2)
  c(mean(zero_term), mean(first_term), mean(second_term), mean(third_term), var_est)
}

write.csv(sim_100_withforest, "~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-15/var_explore_2_1_nof_rf_large_rho7.csv" )


rho = .5

sim_100_withforest = foreach(s = 1:160, .combine = "rbind")%dopar%{
  # TT1 = Sys.time()
  set.seed(s+2023)
  
  n    = 2000
  rho  = 0.5
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X = cbind(X_12, X_310)
  colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y = c(X%*%beta + rnorm(n, 0, 1))
  
  reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  full_model    = randomForest(X[,], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  
  set.seed(s+347)
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X_2 = cbind(X_12, X_310)
  colnames(X_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y_2 = c(X_2%*%beta + rnorm(n, 0, 1))
  
  
  N = 200
  zero_term = c()
  first_term = c()
  second_term = c()
  third_term = c()
  for (i in 1:n){
    beta_curr = beta[-1]
    beta_curr[1] = (1+rho)*beta_curr[1]
    curr_y = rnorm(N, beta_curr%*%X_2[i,-1],sqrt(2-rho^2))
    curr_x = rnorm(N, X_2[i,2]*rho, sqrt(1-rho^2))
    
    second_data =  as.data.frame(cbind(curr_x, matrix(rep(X_2[i,-1], N), ncol = 9, byrow = T)))
    colnames(second_data) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
    
    zero_term[i] = 2*mean((curr_y - predict(full_model, X_2[i,]))*
                            (Y_2[i] - predict(full_model, X_2[i,])))
    first_term[i] = mean((curr_y  - predict(full_model, X_2[i,]))^2)
    second_term[i] = mean((curr_y -  predict(full_model, newdata = second_data))^2)
    third_term[i] = mean((Y_2[i] -  predict(full_model, newdata = second_data))^2)
  }
  #### Variance and Std
  base_perm = mean( first_term - second_term + third_term)
  var_est = mean(( first_term - second_term + third_term-base_perm)^2)
  c(mean(zero_term), mean(first_term), mean(second_term), mean(third_term), var_est)
}

write.csv(sim_100_withforest, "~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-15/var_explore_2_1_nof_rf_large_rho5.csv" )


rho = .3

sim_100_withforest = foreach(s = 1:160, .combine = "rbind")%dopar%{
  # TT1 = Sys.time()
  set.seed(s+2023)
  
  n    = 2000
  rho  = 0.3
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X = cbind(X_12, X_310)
  colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y = c(X%*%beta + rnorm(n, 0, 1))
  
  reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  full_model    = randomForest(X[,], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  
  set.seed(s+347)
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X_2 = cbind(X_12, X_310)
  colnames(X_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y_2 = c(X_2%*%beta + rnorm(n, 0, 1))
  
  
  N = 200
  zero_term = c()
  first_term = c()
  second_term = c()
  third_term = c()
  for (i in 1:n){
    beta_curr = beta[-1]
    beta_curr[1] = (1+rho)*beta_curr[1]
    curr_y = rnorm(N, beta_curr%*%X_2[i,-1],sqrt(2-rho^2))
    curr_x = rnorm(N, X_2[i,2]*rho, sqrt(1-rho^2))
    
    second_data =  as.data.frame(cbind(curr_x, matrix(rep(X_2[i,-1], N), ncol = 9, byrow = T)))
    colnames(second_data) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
    
    zero_term[i] = 2*mean((curr_y - predict(full_model, X_2[i,]))*
                            (Y_2[i] - predict(full_model, X_2[i,])))
    first_term[i] = mean((curr_y  - predict(full_model, X_2[i,]))^2)
    second_term[i] = mean((curr_y -  predict(full_model, newdata = second_data))^2)
    third_term[i] = mean((Y_2[i] -  predict(full_model, newdata = second_data))^2)
  }
  #### Variance and Std
  base_perm = mean( first_term - second_term + third_term)
  var_est = mean(( first_term - second_term + third_term-base_perm)^2)
  c(mean(zero_term), mean(first_term), mean(second_term), mean(third_term), var_est)
}

write.csv(sim_100_withforest, "~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-15/var_explore_2_1_nof_rf_large_rho3.csv" )






rho = .1

sim_100_withforest = foreach(s = 1:160, .combine = "rbind")%dopar%{
  # TT1 = Sys.time()
  set.seed(s+2023)
  
  n    = 2000
  rho  = 0.1
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X = cbind(X_12, X_310)
  colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y = c(X%*%beta + rnorm(n, 0, 1))
  
  reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  full_model    = randomForest(X[,], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  
  set.seed(s+347)
  X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310 = matrix(runif(n*8, 0,1), ncol = 8)
  X_2 = cbind(X_12, X_310)
  colnames(X_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y_2 = c(X_2%*%beta + rnorm(n, 0, 1))
  
  
  N = 200
  zero_term = c()
  first_term = c()
  second_term = c()
  third_term = c()
  for (i in 1:n){
    beta_curr = beta[-1]
    beta_curr[1] = (1+rho)*beta_curr[1]
    curr_y = rnorm(N, beta_curr%*%X_2[i,-1],sqrt(2-rho^2))
    curr_x = rnorm(N, X_2[i,2]*rho, sqrt(1-rho^2))
    
    second_data =  as.data.frame(cbind(curr_x, matrix(rep(X_2[i,-1], N), ncol = 9, byrow = T)))
    colnames(second_data) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
    
    zero_term[i] = 2*mean((curr_y - predict(full_model, X_2[i,]))*
                            (Y_2[i] - predict(full_model, X_2[i,])))
    first_term[i] = mean((curr_y  - predict(full_model, X_2[i,]))^2)
    second_term[i] = mean((curr_y -  predict(full_model, newdata = second_data))^2)
    third_term[i] = mean((Y_2[i] -  predict(full_model, newdata = second_data))^2)
  }
  #### Variance and Std
  base_perm = mean( first_term - second_term + third_term)
  var_est = mean(( first_term - second_term + third_term-base_perm)^2)
  c(mean(zero_term), mean(first_term), mean(second_term), mean(third_term), var_est)
}

write.csv(sim_100_withforest, "~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-15/var_explore_2_1_nof_rf_large_rho1.csv" )
















