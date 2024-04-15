### Trying to replicate the work from Vanelli and Wasserman and see if the TMLE 
### would improve their performance.

library(grf)
library(mgcv)
library(MASS)
library(car)
library(ks)
library(akima)
library(KernSmooth)


## Estimation of the decorrelated LOCO with iteration
## The multivariate estimation can be as slow as it gets, which would be time 
## consuming as well. Thus, I decided to consider both w and z to be 1d to 
## avoid this issue.

set.seed(2013)
rho = 0.3

w = rnorm(1000, 3, 1)
y = 2 *w + rnorm(length(w), 0, 1)
z = rho * w + rnorm(length(w), 0, 1)
x = cbind(w,z)

full_model = lm(y~1+w+z)
reduced_model = lm(y~1+z)

full_predict = cbind(rep(1,dim(x)[1]),x)%*%full_model$coefficients
reduced_predict = cbind(rep(1,length(z)), z)%*%reduced_model$coefficients
Regular_LOCO_hat_1 = mean((full_predict-reduced_predict)^2)

## One-step decorrelated LOCO
## Notice that this whole chunk of code is only used to estimate the 1d
## case
pw = density(w)
pz = density(z)
pwz = kde2d(w, z, n=100)


## We need to be extra cautious here, the prediction function is sensitive
## to the covariate names.

# Returns a numerical value of 
mu_hat_star_z <- function(inputz, n = 100, pw, full_model){
  ## Obtaining samples from p(w)
  ## set.seed(2014)
  samples = sample(pw$x, size = n, prob = pw$y)
  predict_data = data.frame(x = cbind(samples, rep(inputz, length(samples))))
  colnames(predict_data) = c("w", "z")
  result  = mean(predict(full_model, predict_data))
  return(result)
}

mu_hat_star_x <- function(inputx, n = 100, pz, full_model){
  samples = sample(pz$x, size = n, prob = pz$y)
  predict_data = data.frame(x = cbind(rep(inputx, length(samples)), samples))
  colnames(predict_data) = c("w", "z")
  result = mean(predict(full_model, predict_data))
  return(result)
}

second_term_x <- function(inputx, n = 100, pz, full_model, reduced_model){
  #set.seed(2024)
  samples = sample(pz$x, size = n, prob = pz$y)
  predict_data_full = data.frame(x = cbind(rep(inputx, length(samples)), samples))
  colnames(predict_data_full) = c("w", "z")
  predict_data_reduced = data.frame(z = samples)
  #print(head(predict(full_model, predict_data_full)))
  #print(head(predict(reduced_model, predict_data_reduced)))
  result = mean((predict(full_model, predict_data_full)- predict(reduced_model, predict_data_reduced))^2)
  return(result)
}

first_term_z <- function(inputz, n = 100, pw, full_model, reduced_model){
  ## Obtaining samples from p(w)
  ## set.seed(2014)
  samples = sample(pw$x, size = n, prob = pw$y)
  predict_data_full = data.frame(x = cbind(samples, rep(inputz, length(samples))))
  colnames(predict_data_full) = c("w", "z")
  predict_data_reduced = data.frame(z = rep(inputz, length(samples)))
  result = mean((predict(full_model, predict_data_full)- predict(reduced_model, predict_data_reduced))^2)
  return(result)
}

k2d_approx <- function(xout, yout, k2ddensity){
  density_values <- k2ddensity$z
  grid_points <- expand.grid(x = k2ddensity$x, y = k2ddensity$y)
  approximated_density <- suppressWarnings(interp(grid_points$x, grid_points$y, density_values, xout, yout))
  #approximated_density <- suppressWarnings(interp(k2ddensity$x, k2ddensity$y, density_values, xout, yout))
  
  return(approximated_density)
}

density_ratio <- function(inputx, inputz, pw, pz, pwz){
  ratio = approx(pw$x, pw$y, xout = inputx)$y*
  approx(pz$x, pz$y, xout = inputz)$y/
  k2d_approx(inputx, inputz, pwz)$z
  return(ratio)
}
third_term_xz <- function(inputx, inputz, inputy, pw, pz, pwz, full_model, reduced_model){
  result = density_ratio(inputx, inputz, pw, pz, pwz)*
    (c(predict(full_model,data.frame(w=inputx, z = inputz)-mu_hat_star_z(inputx, 100, pz, full_model))))*
    (inputy - predict(full_model,data.frame(w=inputx, z = inputz)))
  return(result)
}

# We (actcually don't) need to vectorize this function if we were to use it in the future, it is too slow
# Mainly because the computation in finding the pwz
# Idea: we can try to improve this function by sampling a list of (x,y)

DecorLOCO <- function(y, w, z, pw, pz, pwz, full_model, reduced_model){
  first_term = c()
  second_term = c()
  third_term = c()
  #time1 = Sys.time()
  for(i in 1:length(y)){
    first_term[i] = first_term_z(z[i], 100, pw, full_model, reduced_model)
    second_term[i] = second_term_x(w[i], 100, pz, full_model, reduced_model)
    third_term[i] = c(third_term_xz(w[i], z[i], y[i], pw, pz, pwz,full_model, reduced_model))
  }
  #time2 = Sys.time()
  #time2 - time1
  return(mean(first_term)/2 + mean(second_term)/2 + mean(third_term))
  
}

# The speed of computing is too slow now.
DecorLOCO(y, w, z, pw, pz, pwz, full_model, reduced_model)

third_term = c()
for(i in 1:length(y)){
       third_term[i] = c(third_term_xz(w[i], z[i], y[i], pw, pz, pwz,full_model, reduced_model))
}


### Testing Region

set.seed(2014)
sample(pw$x, size = 10, prob = pw$y)
set.seed(2014)
sample(px$x, size = 10, prob = px$y)
mu_hat_z(4, 100, pw)



###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


## Generating data
set.seed(2024325)
rho = 0.99
Sig = diag(rep(1,10))
Sig[1,2] = rho
Sig[2,1] = rho
x = mvrnorm(1000, rep(0,10), Sig)
w = x[, 1]
z = x[, -1]


## Simulation process
beta = 5
y = beta*w + rnorm(length(w), 0, 1)


## Create a data frame with y and the previously defined x
colnames(x) = c("w1", "z1", "z2", "z3", "z4","z5","z6","z7","z8","z9")
data_df <- data.frame(cbind(y, x))

## Fit the linear regression model
lmmodel <- lm(formula = y ~1+ w1 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9, data = data_df)
vif(lmmodel)

## Estimate the traditional LOCO
full_model = lm(y~x)
reduced_model = lm(y~z)

full_predict = cbind(rep(1,dim(x)[1]),x)%*%full_model$coefficients
reduced_predict = cbind(rep(1,dim(z)[1]), z)%*%reduced_model$coefficients
Regular_LOCO_hat_1 = mean((full_predict-reduced_predict)^2)
Regular_LOCO_hat_2 =var(y)*(summary(full_model)$r.squared - summary(reduced_model)$r.squared)
Regular_LOCO_hat_2

## Estimation of the decorrelated LOCO
px = density(w)
pz = kde(z)

