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


full_model_lm = lm(Y~0+X)
reduced_model_lm = lm(Y~0+X[,-1])



n = length(Y)
N = 300

# We can now calculate the eif function
zero_term = c()
first_term = c()
second_term = c()
third_term = c()

for(i in 1:n){
  beta_curr = beta[-1]
  beta_curr[1] = (1+rho)*beta_curr[1]
  curr_y =  rnorm(N, beta_curr%*%X[i,-1],2-rho^2)
  
  ### Zero Term
  # zero_term[i] = 2*mean(rnorm(100,0,2-rho^2)*
  # zero_term[i] = zero_term[i] = 2*mean(curr_y*
  zero_term[i] = 2*mean((curr_y- predict(full_model, X[i, ]))*
                                         (Y[i] - predict(full_model, X[i, ])))
  # zero_term[i] =  mean((curr_y- predict(full_model, X[i, ])))
  # zero_term[i] =  mean((curr_y- c(full_model_lm$coefficients%*%X[i,])))
}
mean(zero_term)


set.seed(2023)

n    = 2000
rho  = 0.9
X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
# X_12 = rbinormcop(n, rho)
X_310 = matrix(runif(n*8, 0,1), ncol = 8)
X = cbind(X_12, X_310)
colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)

Y = c(X%*%beta + rnorm(n, 0, 2))


### We'll start with the cpi of the first feature. 
reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                             ntree = 500, keep.inbag = T)
full_model    = randomForest(X[,], Y, nodesize = 5, 
                             ntree = 500, keep.inbag = T)
reduced_den = empirical_rf_pdf(reduced_model, X[,-1], X[,-1])


full_model_lm = lm(Y~0+X)
reduced_model_lm = lm(Y~0+X[,-1])



n = length(Y)
N = 300

# We can now calculate the eif function
zero_term = c()
first_term = c()
second_term = c()
third_term = c()

for(i in 1:n){
  beta_curr = beta[-1]
  beta_curr[1] = (1+rho)*beta_curr[1]
  curr_y =  rnorm(N, beta_curr%*%X[i,-1],2-rho^2)
  
  ### Zero Term
  # zero_term[i] = 2*mean(rnorm(100,0,2-rho^2)*
  # zero_term[i] = zero_term[i] = 2*mean(curr_y*
  zero_term[i] = 2*mean((curr_y- predict(full_model, X[i, ]))*
                          (Y[i] - predict(full_model, X[i, ])))
  # zero_term[i] =  mean((curr_y- predict(full_model, X[i, ])))
  # zero_term[i] =  mean((curr_y- c(full_model_lm$coefficients%*%X[i,])))
}
mean(zero_term)

set.seed(2023)

n    = 2000
rho  = 0.9
X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
# X_12 = rbinormcop(n, rho)
X_310 = matrix(runif(n*8, 0,1), ncol = 8)
X = cbind(X_12, X_310)
colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)

Y = c(X%*%beta + rnorm(n, 0, 1/2))


### We'll start with the cpi of the first feature. 
reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                             ntree = 500, keep.inbag = T)
full_model    = randomForest(X[,], Y, nodesize = 5, 
                             ntree = 500, keep.inbag = T)
reduced_den = empirical_rf_pdf(reduced_model, X[,-1], X[,-1])


full_model_lm = lm(Y~0+X)
reduced_model_lm = lm(Y~0+X[,-1])



n = length(Y)
N = 300

# We can now calculate the eif function
zero_term = c()
first_term = c()
second_term = c()
third_term = c()

for(i in 1:n){
  beta_curr = beta[-1]
  beta_curr[1] = (1+rho)*beta_curr[1]
  curr_y =  rnorm(N, beta_curr%*%X[i,-1],2-rho^2)
  
  ### Zero Term
  # zero_term[i] = 2*mean(rnorm(100,0,2-rho^2)*
  # zero_term[i] = zero_term[i] = 2*mean(curr_y*
  zero_term[i] = 2*mean((curr_y- predict(full_model, X[i, ]))*
                          (Y[i] - predict(full_model, X[i, ])))
  # zero_term[i] =  mean((curr_y- predict(full_model, X[i, ])))
  # zero_term[i] =  mean((curr_y- c(full_model_lm$coefficients%*%X[i,])))
}
mean(zero_term)

a = rnorm(300, 0, 1)
b = rnorm(300, 0,1)
mean(a*b)

