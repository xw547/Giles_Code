library(MASS)
library(mvtnorm)
library(doParallel)
library(ks)
library(ggplot2)
library(tidyr)
registerDoParallel(8)

### We may start with a simple logistic regression model

n = 1000
set.seed(2023)

Sigma = .5^(2-toeplitz(2:1))
X = mvrnorm(n, mu = c(0, 0), Sigma)
prop = 1- 1/(1 + exp(- X[,1] + .5* X[,2]))
treat = rbinom(n, 1, prop)

model_hat = glm(treat~X, family = "binomial")
prop_hat = model_hat$fitted.values
bias = abs(1/prop - 1/prop_hat)





### Now we'll proceed to the density ratio model

w = X[, 1]
z = X[, 2]
x = cbind(w,z)

### Density estimation:
### Notice that we are only interested in the density estimation of existing
### points.

pw = density(w)
pz = density(z)

### The following kde density estimates were used to estimate the density ratio 
### in the third term.

pw_index = kde(w, eval.points = w, density = T)
pz_index = kde(z, eval.points = z, density = T)
pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)


### The speed of computing can be greatly accelerated with prespecified eval.
### points.

### We can replace index with the which if needed.

density_ratio <- function(index, pw_index, pz_index, pwz){
  curr_pw  = pw_index$estimate[index]
  curr_pz  = pz_index$estimate[index]
  curr_pwz = pwz$estimate[index]
  return(curr_pw*curr_pz/curr_pwz)
}


oracle_pw_density <- function(x){
  dnorm(x, 0, 1)
}

oracle_pz_density <- function(x){
  dnorm(x, 0, 1)
}

oracle_pwz_density <- function(input){
  sigma_matrix <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
  mean_vector <- c(0, 0)
  dmvnorm(input, mean = mean_vector, sigma = sigma_matrix)
}

oracle_density_ratio <- function(w,z){
  result = oracle_pw_density(w)*oracle_pz_density(z)/oracle_pwz_density(c(w,z))
  return(result)
}


test = foreach(i = 1:length(w), .combine = "c")%dopar%{
  abs(density_ratio(i, pw_index, pz_index, pwz) - oracle_density_ratio(w[i],z[i]) )
}

quantile(bias)
quantile(test)
