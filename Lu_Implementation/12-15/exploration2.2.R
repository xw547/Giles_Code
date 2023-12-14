### The idea would be to use sample splitting 
### technique to reduce the bias. 

### Loading required packages
library("MASS")
library("randomForest")
library("MASS")
source("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/11-02/empirical_weight.R")
library(doParallel)
library(foreach)
registerDoParallel(7)

set.seed(2023)

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
# Speeeeeed reduced_den = empirical_rf_pdf(reduced_model, X[,-1], X[,-1])


set.seed(2023+347)
X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
# X_12 = rbinormcop(n, rho)
X_310 = matrix(runif(n*8, 0,1), ncol = 8)
X_2 = cbind(X_12, X_310)
colnames(X_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)

Y_2 = c(X_2%*%beta + rnorm(n, 0, 1))

reduced_model_2 = randomForest(X_2[,-1], Y_2, nodesize = 5, 
                             ntree = 500, keep.inbag = T)
full_model_2    = randomForest(X_2[,], Y_2, nodesize = 5, 
                             ntree = 500, keep.inbag = T)


N = 200
zero_term = c()
first_term = c()
second_term = c()
third_term = c()
zero_term_2 = c()
first_term_2 = c()
second_term_2 = c()
third_term_2 = c()
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
  
  curr_y_2 = rnorm(N, beta_curr%*%X[i,-1],sqrt(2-rho^2))
  curr_x_2 = rnorm(N, X[i,2]*rho, sqrt(1-rho^2))
  second_data_2 =  as.data.frame(cbind(curr_x_2, matrix(rep(X[i,-1], N), ncol = 9, byrow = T)))
  colnames(second_data_2) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
  
  zero_term_2[i] = 2*mean((curr_y_2 - predict(full_model_2, X[i,]))*
                          (Y[i] - predict(full_model_2, X[i,])))
  first_term_2[i] = mean((curr_y_2  - predict(full_model_2, X[i,]))^2)
  second_term_2[i] = mean((curr_y_2 -  predict(full_model_2, newdata = second_data_2))^2)
  third_term_2[i] = mean((Y[i] -  predict(full_model_2, newdata = second_data_2))^2)
}

mean(zero_term)
mean(first_term)
mean(second_term)
mean(third_term)
mean(zero_term_2)
mean(first_term_2)
mean(second_term_2)
mean(third_term_2)
mean(c(zero_term, zero_term_2))
mean(c(first_term, first_term_2))
mean(c(second_term, second_term_2))
mean(c(third_term, third_term_2))

