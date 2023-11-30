### This is a package used to calculate the density ratio, 
### which used to be the speed constraint in the 
### EIF estimation.

### library(KernSmooth): The eval.points can't be determined.
library(fields)
library(ks)

### Data Preparation

set.seed(2013)
rho = 0.3

w = rnorm(1000, 3, 1)
y = 2 *w + rnorm(length(w), 0, 1)
z = rho * w + rnorm(length(w), 0, 1)
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

mu_hat_star_z <- function(inputz, n = 100, pw, full_model){
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









