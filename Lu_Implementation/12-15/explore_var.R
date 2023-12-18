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
ve = rnorm(n)

set.seed(20394)
X_12 = mvrnorm(n, rep(0.5,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
# X_12 = rbinormcop(n, rho)
X_310 = matrix(runif(n*8, 0,1), ncol = 8)
X_2 = cbind(X_12, X_310)


err = ve%*%(X)%*%solve(t(X)%*%(X))%*%X[i,]%*%t(X[i,])%*%solve(t(X)%*%(X))%*%t(X)%*%ve
err_2 = ve%*%(X)%*%solve(t(X)%*%(X))%*%X_2[i,]%*%t(X_2[i,])%*%solve(t(X)%*%(X))%*%t(X)%*%ve
err_2_full = ve%*%(X)%*%solve(t(X)%*%(X))%*%t(X_2)%*%(X_2)%*%solve(t(X)%*%(X))%*%t(X)%*%ve
