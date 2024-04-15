library(MASS)
n = 2000000
sig = diag(rep(1,10))


nonlinear = foreach(rho = seq(.1, .9, .1), .combine = "c")%do%{
  sig[1,2] = rho
  sig[2,1] = rho
  
  X = rmvnorm(n, mean = rep(0,10), sigma = sig)
  mean((lm(sin(X[,1])~X[,2])$residuals)^2)*25
}

linear = 25*(1-seq(.1, .9, .1)^2)
nonlinear[1] = (1 - exp(-2))/2*25

nonlinear_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-03-11/nonlinear.csv")[, -1]
linear_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-03-11/linear.csv")[, -1]


estimates_linear = matrix(rowMeans(linear_data[sort(c(2+4*(0:8), 3+4*(0:8))),]), byrow = F, nrow = 2)
estimates_nonlinear = matrix(rowMeans(nonlinear_data[sort(c(2+4*(0:8), 3+4*(0:8))),]), byrow = F, nrow = 2)



