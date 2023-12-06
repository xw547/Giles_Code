### In this document, I'll specify the 

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

full_model_lm = lm(Y~0+X)
reduced_model_lm = lm(Y~0+X[,-1])

N = 200
zero_term = c()
first_term = c()
second_term = c()
third_term = c()
for (i in 1:n){
  beta_curr = beta[-1]
  beta_curr[1] = (1+rho)*beta_curr[1]
  curr_y = rnorm(N, beta_curr%*%X[i,-1],sqrt(2-rho^2))
  curr_x_1 = rnorm(N, X[i,2]*rho, sqrt(1-rho^2))
  curr_x  = cbind(curr_x_1, matrix(X[i, -1], nrow = N, ncol = 9, byrow = T))
  
  zero_term[i] = 2*mean((curr_y - full_model_lm$fitted.values[i])*
                          (Y[i] - full_model_lm$fitted.values[i]))
  first_term[i] = mean((curr_y  - full_model_lm$fitted.values[i])^2)
  second_term[i] = mean((curr_y -  curr_x%*%full_model_lm$coefficients )^2)
  third_term[i] = mean((Y[i] -  curr_x%*%full_model_lm$coefficients )^2)
}
mean(zero_term)
mean(first_term)
mean(second_term)
mean(third_term)

