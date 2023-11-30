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

set.seed(2025)

n    = 2000
rho  = 0.5
X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
# X_12 = rbinormcop(n, rho)
X_310 = matrix(runif(n*8, 0,1), ncol = 8)
X = cbind(X_12, X_310)
colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)

Y = c(X%*%beta + rnorm(n, 0, .1))
full_model_lm = lm(Y~0+X)
full_model = randomForest(X, Y)

N = 200
zero_term = c()
first_term = c()

for(i in 1:n){
  first_term[i] = t(beta) %*% (X[i,]%*%t(X[i,]) - colMeans(X)%*% t(colMeans(X))) %*% beta
}

for(i in 1:n){
  curr_y = c(beta%*%X[i,]) + rnorm(N)
  beta_curr = beta[-1]
  beta_curr[1] = (1+rho)*beta_curr[1]
  curr_y_2 =  rnorm(N, beta_curr%*%X[i,-1],2-rho^2)
  zero_term[i] = mean((Y[i]  - full_model$predicted[i]) * (Y[i] - full_model$predicted[i]))
  first_term[i] = mean(curr_y_2 * (Y[i] - full_model$predicted[i]))
  #first_term[i] = mean(curr_y_2/curr_y)
  #zero_term[i] = c(beta%*%X[i,]) * (Y[i]- full_model$predicted[i])
  #zero_term[i] = c(beta%*%X[i,]) * (Y[i]- predict(full_model, X[i, ]))
  #zero_term[i] = c(beta%*%X[i,]) * c(beta%*%(X[i,] - colMeans(X[-i,])))
}
mean(zero_term)
mean(first_term)

for(i in 1:n){
  beta_curr = beta[-1]
  beta_curr[1] = (1+rho)*beta_curr[1]
  curr_y =  rnorm(N, beta_curr%*%X[i,-1],1)
  
  ### Zero Term
  zero_term[i] = 2*mean(curr_y -full_model_lm$fitted.values[i])*
    (Y[i] - full_model_lm$fitted.values[i])
  #first_term[i] = mean(curr_y - full_model$predicted)
}
mean(zero_term)


