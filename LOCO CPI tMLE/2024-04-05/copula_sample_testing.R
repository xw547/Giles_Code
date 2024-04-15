rmvcopula <- function(n, p, Sigma){
  require("mvtnorm")
  A = chol(Sigma)
  samples = rmvnorm(n, rep(0, p), diag(rep(1, p)))
  return(pnorm(samples %*% A))
}
set.seed(2024)


copula_ff = c()
for (i in 1:9){
  rho =  0.1*i

  
  sig = diag(rep(1,10))
  sig[1,2] = rho
  sig[2,1] = rho
  n_1 = 1e6
  
  X = rmvcopula(n_1, 10, sig)
  X_reduced = X[,-1]
  beta = c(5, rep(0,9))
  Y = X%*%beta 
  
  
  full_data = data.frame(y = Y, X = X)
  colnames(full_data) <- c("output", paste0("X",1:10))
  reduced_data = data.frame(y = Y, X_reduced = X_reduced)
  colnames(reduced_data) <- c("output", paste0("X",2:10))
  
  
  X1.X.2_model = lm(X1~X2, data = full_data[1:n_1,])

  copula_ff[i] = 25*mean((full_data[1:n_1,2] -X1.X.2_model$fitted.values)^2)
  
}
