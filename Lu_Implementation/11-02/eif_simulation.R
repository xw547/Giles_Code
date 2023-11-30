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

sim_100_withforest = foreach(s = 1:160, .combine = "cbind")%dopar%{
  set.seed(s+2023)
  # TT1 = Sys.time()
  n    = 2000
  rho  = 0.5
  X_12 = mvrnorm(n, rep(0,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
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
  
  
  
  n = length(Y)
  N = 100
  
  # We can now calculate the eif function
  zero_term = c()
  first_term = c()
  second_term = c()
  third_term = c()
  
  for(i in 1:n){
    curr_y = Y[sample(1:n, N, replace = T, prob = reduced_den[i,])]
    curr_x = X[sample(1:n, N, replace = T, prob = reduced_den[i,]),1]
    
    ### Zero Term
    zero_term = 2*(predict(reduced_model, X[i, -1]) - predict(full_model, X[i, ]))*
      (Y[i] - predict(full_model, X[i, ]))
    
    ### First Term
    curr_pred = predict(full_model, X[i, ])
    first_term[i] = mean((curr_y-curr_pred)^2)
    
    ### Second Term
    second_data =  as.data.frame(cbind(curr_x, matrix(rep(X[i,-1], N), ncol = 9, byrow = T)))
    colnames(second_data) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
    second_term[i] = mean((curr_y%*%t(rep(1,N))- t(predict(full_model, second_data)%*%t(rep(1,N))))^2)
    
    ### Third Term
    third_pred = predict(full_model, second_data)
    third_term[i] = mean((Y[i] - predict(full_model, second_data))^2)
  }
  
  eif_estimate = mean(zero_term + first_term - second_term + third_term)
  empirical_eif = 2*mean((lm(X[,1]~X[,2])$residuals)^2)
  
  c(eif_estimate, empirical_eif)
}
sim_100_withforest = read.csv("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/sim_160_rho_5.csv", sep = " ")

eif_results = as.numeric(sim_100_withforest[1,])
t_95 = qt(.975, 160-1)
eif_mean = 1.5
SE = sqrt(var(eif_results))
eif_Lower = eif_mean - t_95*SE
eif_Upper = eif_mean + t_95*SE



