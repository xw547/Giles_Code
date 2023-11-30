library(mvtnorm)
library(doParallel)
library(ks)
library(ggplot2)
library(tidyr)
registerDoParallel(8)


set.seed(2013)
rho = 0

n = 10000

w = rnorm(n, 3, 1)
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


oracle_pw_density <- function(x){
  dnorm(x, 3, 1)
}

oracle_pz_density <- function(x){
  dnorm(x, 0, 1)
}

oracle_pwz_density <- function(input){
  sigma_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
  mean_vector <- c(3,0)
  dmvnorm(input, mean = mean_vector, sigma = sigma_matrix)
}

oracle_density_ratio <- function(w,z){
  result = oracle_pw_density(w)*oracle_pz_density(z)/oracle_pwz_density(c(w,z))
  return(result)
}


test = foreach(i = 1:length(w), .combine = "c")%dopar%{
  abs(density_ratio(i, pw_index, pz_index, pwz) - oracle_density_ratio(w[i],z[i]) )
}

length(which(test>10/sqrt(length(w))))

test_two = foreach(i = 1:length(w), .combine = "c")%dopar%{
  abs(pw_index$estimate[i] - oracle_pw_density(w[i]) )
}

length(which(test_two>10/sqrt(length(w))))

test_three = foreach(i = 1:length(w), .combine = "c")%dopar%{
  abs(pwz$estimate[i]- oracle_pwz_density(c(w[i],z[i])))
}

ggplot(data = data.frame(Value = test), aes(x = 1, y = Value)) +
  geom_boxplot() +
  labs(x = "", y = "Test Values") +
  ggtitle("Boxplot of Test Values")


test_df <- as.data.frame(cbind(test, test_two, test_three))
colnames(test_df) <- c("Ratio", "1d", "2d" )  # Replace with your column names

# Reshape data into long format
test_long <- gather(test_df, key = "Column", value = "Value")

# Create a facetted boxplot using ggplot
ggplot(test_long, aes(x = Column, y = Value)) +
  geom_boxplot() +
  labs(x = "Column", y = "Values") +
  ggtitle("Boxplots of absolute value of bias") +
  facet_wrap(~ Column, scales = "free")

