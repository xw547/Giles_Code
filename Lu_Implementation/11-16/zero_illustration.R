### Loading required packages
library("MASS")
library("VGAM")
library("permimp")
library("party")
library("randomForest")
library("MASS")
source("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/11-02/empirical_weight.R")
library(doParallel)
library(foreach)
registerDoParallel(7)


zero_collection <- foreach(s = 1: 80, .combine = "c")%dopar%{
  
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
  
  ### We'll start with the cpi of the first feature. 
  reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  full_model    = randomForest(X[,], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  reduced_den = empirical_rf_pdf(reduced_model, X[,-1], X[,-1])
  
  full_model = randomForest(X, Y)
  reduced_model = randomForest(X[,-1], Y)
  
  
  full_model_lm = lm(Y~X)
  reduced_model_lm = lm(Y~X[,-1])
  
  
  
  n = length(Y)
  N = 200
  
  # We can now calculate the eif function
  zero_term = c()
  for(i in 1:n){
    beta_curr = beta[-1]
    beta_curr[1] = (1+rho)*beta_curr[1]
    curr_y =  rnorm(N, beta_curr%*%X[i,-1],2-rho^2)
    
    ### Zero Term
    zero_term[i] = 2*mean(curr_y - predict(full_model, X[i, ]))*
      (Y[i] - predict(full_model, X[i, ]))
    #first_term[i] = mean(curr_y - full_model$predicted)
  }
  mean(zero_term)
}


zero_collection_xfitting <- foreach(s = 1: 80, .combine = "c")%dopar%{
  
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
  
  ### We'll start with the cpi of the first feature. 
  reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  full_model    = randomForest(X[,], Y, nodesize = 5, 
                               ntree = 500, keep.inbag = T)
  reduced_den = empirical_rf_pdf(reduced_model, X[,-1], X[,-1])
  
  full_model = randomForest(X, Y)
  reduced_model = randomForest(X[,-1], Y)
  
  
  full_model_lm = lm(Y~X)
  reduced_model_lm = lm(Y~X[,-1])
  
  n    = 2000
  rho  = 0.9
  X_12_2 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
  # X_12 = rbinormcop(n, rho)
  X_310_2 = matrix(runif(n*8, 0,1), ncol = 8)
  X_2= cbind(X_12_2, X_310_2)
  colnames(X_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)
  
  Y_2 = c(X_2%*%beta + rnorm(n, 0, 1))
  
  n = length(Y_2)
  N = 200
  
  # We can now calculate the eif function
  zero_term = c()
  for(i in 1:n){
    beta_curr = beta[-1]
    beta_curr[1] = (1+rho)*beta_curr[1]
    curr_y =  rnorm(N, beta_curr%*%X_2[i,-1],2-rho^2)
    
    ### Zero Term
    zero_term[i] = 2*mean(curr_y - predict(full_model, X_2[i, ]))*
      (Y_2[i] - predict(full_model, X_2[i, ]))
    #first_term[i] = mean(curr_y - full_model$predicted)
  }
  mean(zero_term)
}