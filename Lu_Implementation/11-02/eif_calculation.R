### This is the function to calculate the EIF of CPI with a empirical dist
### For now, the empirical dist was obtained through the LuHardin Paper
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

set.seed(2022)
n    = 2000
rho  = 0.5
X_12 = mvrnorm(n, rep(0,2), Sigma = matrix(c(1, rho, rho,1), ncol =2))
#X_12 = rbinormcop(n, rho)

X_310 = matrix(runif(n*8, 0,1), ncol = 8)
X = cbind(X_12, X_310)
colnames(X) <- c("1", "2", "3", "4", "5", "6", "7","8", "9", "10")
beta = c(1, 1, 1, 1, 1, 0, 0.5, 0.8, 1.2, 1.5)

Y = c(X%*%beta + rnorm(n, 0, 1))

### We'll start with the required models
reduced_model = randomForest(X[,-1], Y, nodesize = 5, 
                             ntree = 500, keep.inbag = T)
full_model    = randomForest(X[,], Y, nodesize = 5, 
                             ntree = 500, keep.inbag = T)
reduced_den = empirical_rf_pdf(reduced_model, X[,-1], X[,-1])

### We can now consider the difference in them. 
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



