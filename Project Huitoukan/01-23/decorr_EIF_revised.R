### EIF calculation 
### In this R file, I'll try to reproduce the result from Wasserman first. 
### This is a package used to calculate the density ratio, 
### which used to be the speed constraint in the 
### EIF estimation.

### library(KernSmooth): The eval.points can't be determined.
library(fields)
library(ks)

### Data generation

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

pz_index_2 = kde(z_2, eval.points = z_2, density = T)
pwz_2 = kde(cbind(w_2,z_2), eval.points = cbind(w_2,z_2), density = T)


### In this version, we'll 



