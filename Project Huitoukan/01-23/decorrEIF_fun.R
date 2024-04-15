### EIF calculation 
### In this R file, I'll try to reproduce the result from Wasserman first. 
### This is a package used to calculate the density ratio, 
### which used to be the speed constraint in the 
### EIF estimation.

### library(KernSmooth): The eval.points can't be determined.
library(fields)
library(ks)

### Data generation: In principle, we should be using the Friedman equation 
### to estimate the stuff. 

set.seed(2013)
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
### estimation. However, one thing to mention is that we chose to use index so 
### we may have a better speed. 

pw_index = kde(w, eval.points = w, density = T)
pz_index = kde(z, eval.points = z, density = T)
pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)

### Cross-fitting to apply the normalization.
### Data generation: In principle, we should be using the Friedman equation 
### to estimate the stuff. 

set.seed(2014)
rho = 0.3

w_2 = rnorm(1000, 3, 1)
y_2 = 2 *w_2 + rnorm(length(w_2), 0, 1)
z_2 = rho * w_2 + rnorm(length(w_2), 0, 1)
x_2 = as.data.frame(cbind(w_2,z_2))
colnames(x_2) <- c("w","z")

full_model_2 = lm(y_2~1+w+z,data = x_2)
reduced_model_2 = lm(y_2~1+z, data = x_2)

### Density estimation:
### Notice that we are only interested in the density estimation of existing
### points.

pw_2 = density(w_2)
pz_2 = density(z_2)

### The following kde density estimates were used to estimate the density ratio 
### estimation. However, one thing to mention is that we chose to use index so 
### we may have a better speed. 

pw_index_2 = kde(w_2, eval.points = w_2, density = T)
pw_index_2_1 = kde(w_2, eval.points = w, density = T)
pz_index_2 = kde(z_2, eval.points = z_2, density = T)
pz_index_2_1 = kde(z_2, eval.points = z, density = T)
pwz_2 = kde(cbind(w_2,z_2), eval.points = cbind(w_2,z_2), density = T)
pwz_2_1 = kde(cbind(w_2,z_2), eval.points = cbind(w,z), density = T)



### 

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

density_ratio <- function(index, pw_index, pz_index, pwz){
  curr_pw  = pw_index$estimate[index]
  curr_pz  = pz_index$estimate[index]
  curr_pwz = pwz$estimate[index]
  return(curr_pw*curr_pz/curr_pwz)
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
  result = density_ratio(index, pw_index, pz_index, pwz)*
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


### Now, we proceed to the calculatin of \hat{\psi}
hat_psi <- function(w, z, y, pw, pw_index, 
                    pz, pz_index, pwz, full_model, reduced_model){
  sample_size = length(y)
  phi_vec = c()
  phat_vec = c()
  for (index in 1:sample_size){
    phi_vec[index] <- EIF(index, w[index], z[index], y[index], pw, pw_index, 
                          pz, pz_index, pwz, full_model, reduced_model, psi_hat(n = 100, pw, pz, full_model, reduced_model))
    phat_vec[index] <- pwz$estimate[index]
  }
  phi_vec[abs(phi_vec)>100] = 0
  #print(max(phi_vec))
  #print(min(phi_vec))
  #print(mean(phat_vec*exp(epsilon*phi_vec)*phi_vec))
  #print(mean(phat_vec*exp(epsilon*phi_vec)))
  
  result <-mean(phi_vec) - mean(phat_vec*phi_vec)
  return(result)
}


likelihood_nor <- function(gamma){
  likelihood_individual = c()
  for (index in 1:length(y)){
    likelihood_individual[index] <- pwz$estimate[index]*exp(gamma*tilde_phi[index]) 
    - gamma*tilde_phi[index]
  }
  return(mean(likelihood_individual))
}


### We can start with the revised likelihood function and try to optimize over
### it. 

likelihood_ori <- function(gamma){
  likelihood_individual = c()
  for (index in 1:length(y)){
    likelihood_individual[index] <- pwz$estimate[index]*exp(gamma*eif_vec[index]) 
    - gamma*eif_vec[index]
  }
  return(mean(likelihood_individual))
}