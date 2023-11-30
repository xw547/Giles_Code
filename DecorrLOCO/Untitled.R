### Trying to replicate the work from Vanelli and Wasserman and see if the TL 
### would improve their performance.

library(grf)
library(mgcv)
library(MASS)
library(car)
library(ks)


## Estimation of the decorrelated LOCO with iteration
## The multivariate estimation can be as slow as it gets, which would be time 
## consuming as well. Thus, I decided to consider both w and z to be 1d to 
## avoid this. 

set.seed(2014)
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
pwz = kde2d(w, z, n=50)


## We need to be extra cautious here, the prediction function is sensitive
## to the covariate names.
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
  colnames(predict_data) = c("w", "z")
  predict_data_reduced = data.frame(rep(inputz, length(samples)))
  result = mean((predict(full_model, predict_data_full)- predict(reduced_model, predict_data_reduced))^2)
  return(result)
}


density_ratio = approx(pw$x, pw$y, xout = inputx)*
                approx(pz$x, pz$y, zout = inputz)/
                approx(pwz$)

third_term = 2*mean()
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

