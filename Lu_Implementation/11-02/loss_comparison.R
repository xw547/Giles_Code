`### In this file, I'll try to implement the algorithm mentioned in Lu, 2018 JMLR
### to estimate the joint density of x and y. 

### Loading required packages
library("MASS")
library("VGAM")
library("permimp")
library("party")
library("randomForest")
library("MASS")
source("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/empirical_weight.R")
library(doParallel)
library(foreach)
registerDoParallel(7)
### Problem Set-up
### I'll largely follow the formulation of Hooker and Zhou 2021 SC 

set.seed(2021)
n    = 2000
rho  = 0.5
X_12 = rbinormcop(n, rho)
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

# \int L\left( y, \hat{y}(\tilde{x}, \tilde{z})\right)p(y|\tilde{z})dy
# \sum^n_i\sum^N_j (y_j - \hat(x_i, z_i))^2 
# \int L\left( y, \hat{y}({x}', \tilde{z}) \right)p(y|\tilde{z}) p(x|\tilde{z})dx'dy
# \sum^n_i\sum^N_j\sum^N_t (y_j - \hat(x_j, z_i))^2 
# \int L\left( \tilde{y}, \hat{y}({x}', \tilde{z})\right)p(x'|\tilde{z})dx'
# Notice that we'll be considering the squared error case.

n = length(Y)
N = 100

# We can start with the first term.
first_term = c()
second_term = c()
third_term = c()
TT1 = Sys.time()
for(i in 1:n){
  curr_y = Y[sample(1:n, N, replace = T, prob = reduced_den[i,])]
  curr_x = X[sample(1:n, N, replace = T, prob = reduced_den[i,]),1]
  
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
TT2 = Sys.time()
TT2-TT1

final_first_term = mean(first_term)
final_second_term = mean(second_term)
final_third_term = mean(third_term)
c(final_first_term, final_second_term, final_third_term)





