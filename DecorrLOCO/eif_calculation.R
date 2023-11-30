### This is a package used to calculate the density ratio, 
### which used to be the speed constraint in the 
### EIF estimation.

### library(KernSmooth): The eval.points can't be determined.
library(fields)
library(ks)

### Data generation: In principle, we should be using the Friedman equation 
### to estimate the stuff. 

#set.seed(2013)
rho = 0.3

w = rnorm(1000, 3, 1)
y = 2 *w + rnorm(length(w), 0, 1)
z = rho * w + rnorm(length(w), 0, 1)
x = cbind(w,z)

full_model = lm(y~1+w+z)
reduced_model = lm(y~1+z)

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

## We can replace index with the (which function) if needed.

density_ratio <- function(index, inputw, inputz, pw, pz, pwz){
  curr_pw  = approx(pw$x, pw$y, xout = inputw)$y
  curr_pz  = approx(pz$x, pz$y, xout = inputz)$y
  curr_pwz = pwz$estimate[index]
  return(curr_pw*curr_pz/curr_pwz)
}

## Auxiliary functions for our calculation 
## These functions are inherited from the linear model example, which should've
## been updated.

mu_hat_star_z <- function(inputz, n = 100, pw, full_model){
  samples = sample(pw$x, size = n, prob = pw$y)
  predict_data = data.frame(x = cbind(samples, rep(inputz, length(samples))))
  colnames(predict_data) = c("w", "z")
  result  = mean(predict(full_model, predict_data))
  return(result)
}

mu_hat_star_w <- function(inputw, n = 100, pz, full_model){
  samples = sample(pz$x, size = n, prob = pz$y)
  predict_data = data.frame(x = cbind(rep(inputw, length(samples)), samples))
  colnames(predict_data) = c("w", "z")
  result = mean(predict(full_model, predict_data))
  return(result)
}

first_term_z <- function(inputz, n = 100, pw, full_model, reduced_model){
  samples = sample(pw$x, size = n, prob = pw$y)
  predict_data_full = data.frame(x = cbind(samples, rep(inputz, length(samples))))
  colnames(predict_data_full) = c("w", "z")
  predict_data_reduced = data.frame(z = rep(inputz, length(samples)))
  result = mean((predict(full_model, predict_data_full)-
                   predict(reduced_model, predict_data_reduced))^2)
  return(result)
}

second_term_w <- function(inputw, n = 100, pz, full_model, reduced_model){
  samples = sample(pz$x, size = n, prob = pz$y)
  predict_data_full = data.frame(x = cbind(rep(inputw, length(samples)), samples))
  colnames(predict_data_full) = c("w", "z")
  predict_data_reduced = data.frame(z = samples)
  result = mean((predict(full_model, predict_data_full)-
                   predict(reduced_model, predict_data_reduced))^2)
  return(result)
}

third_term_wz <- function(index, inputw, inputz, inputy, pw, pw_index, 
                          pz, pz_index, pwz, full_model, reduced_model){
  result = density_ratio(index, inputw, inputz, pw, pz, pwz)*
    (c(predict(full_model,data.frame(w=inputw, z = inputz)-
                 mu_hat_star_z(inputw, 100, pz, full_model))))*
    (inputy - predict(full_model,data.frame(w=inputw, z = inputz)))
  return(result)
}

psi_hat <- function(n = 100, pw, pz, full_model, reduced_model){
  sample_w = sample(pw$x, size = n, prob = pw$y)
  sample_z = sample(pz$x, size = n, prob = pz$y)
  predict_data_full = data.frame(x = cbind(sample_w, sample_z))
  colnames(predict_data_full) = c("w", "z")
  predict_data_reduced = data.frame(z = sample_z)
  result = mean((predict(full_model, predict_data_full)-
                   predict(reduced_model, predict_data_reduced))^2)
  return(result)
}

DecorrLOCO <- function(index, inputw, inputz, inputy, pw, pw_index, 
                       pz, pz_index, pwz, full_model, reduced_model){
  third_term = third_term_wz(index, inputw, inputz, inputy, pw, pw_index, 
                            pz, pz_index, pwz, full_model, reduced_model)
  first_term = first_term_z(inputz, n = 100, pw, full_model, reduced_model)
  second_term = second_term_w(inputw, n = 100, pw, full_model, reduced_model)
  return(third_term+ first_term+ second_term)
}

EIF <- function(index, inputw, inputz, inputy, pw, pw_index, 
                pz, pz_index, pwz, full_model, reduced_model, psi_now){
  result = DecorrLOCO(index, inputw, inputz, inputy, pw, pw_index, 
                      pz, pz_index, pwz, full_model, reduced_model) - 2*psi_now
  return(result)
}
### Now, we shall proceed to the computation of EIFs. 
### Requirments:
###             pw, pz, pwz, full_model, reduced_model
###             w, z.

## Making sure there isn't any error.
index = 4
DecorrLOCO(index, w[index], z[index], y[index], pw, pw_index, 
           pz, pz_index, pwz, full_model, reduced_model)

EIF(index, w[index], z[index], y[index], pw, a, 
           pz, b, pwz, full_model, reduced_model, psi_hat(n = 100, pw, pz, full_model, reduced_model))

psi_hat(n = 100, pw, pz, full_model, reduced_model)
### Now, we are able to calculate the EIF efficiently, what remains to be done 
### would be to write necessary code to alternate the density.

### In each iteration, both marginal density will be estimated by 
pwz_grid = kde(cbind(w,z), density = T, gridsize = c(128L, 128L))


## Notice that we have to use the kde2d function to calculate the 
## the joint pdf if we were to do so.
margin_calculate_w <- function(pwz){
  margin_result = list()
  margin_result$x = pwz$eval.points[[1]]
  margin_result$y = rowMeans(pwz$estimate)
  return(margin_result)
 }

margin_calculate_z <- function(pwz){
  margin_result = list()
  margin_result$x = pwz$eval.points[[2]]
  margin_result$y = colMeans(pwz$estimate)
  return(margin_result)
}

likelihood <- function(epsilon, w, z, y, pw, pw_index, 
                pz, pz_index, pwz, pwz_grid, full_model, reduced_model, psi_now){
  margin_w = margin_calculate_w(pwz_grid)
  margin_z = margin_calculate_z(pwz_grid)
  sample_w = sample(x=margin_w$x, prob = margin_w$y, size = length(w), replace = T)
  sample_z = sample(x=margin_z$x, prob = margin_z$y, size = length(w), replace = T)
  eif_vec = c()
  C_epsilon_vec = c()
  psi_hat_now = psi_hat(n = 100, pw, pz, full_model, reduced_model)
  for(index in 1:length(w)){
    curr_eif = EIF(index, w[index], z[index], y[index], pw, a, 
                   pz, b, pwz, full_model, reduced_model, psi_hat_now)
    eif_vec[index] = curr_eif
    C_epsilon_vec[index] = EIF(index, sample_w[index], sample_z[index], y[index], pw, a, 
                        pz, b, pwz, full_model, reduced_model, psi_hat_now)
  }
  result = mean(epsilon*eif_vec) + log(sum(exp(epsilon*C_epsilon_vec)))
  return(result)
}

argmin <- function(epsilon){
  return(likelihood(epsilon, w, z, y, pw, pw_index, 
                                pz, pz_index, pwz, pwz_grid, full_model, reduced_model, psi_now))
}

epsilon_now = optimize(argmin, c(-0.5,0.5))$minimum

if(epsilon_now<1e-5){print("last iteration")} 

augment_margin <- function(pwz, epsilon_now){
  pwz$estimate = 
}
