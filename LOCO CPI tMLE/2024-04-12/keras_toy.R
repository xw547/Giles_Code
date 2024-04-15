library(keras)
library(mvtnorm)
library(ranger)
# Create a synthetic dataset
set.seed(42) # For reproducibility
TT1= Sys.time()
n = 1e3

X_P = cbind(rep(1, n), rmvnorm(n, rep(0, 4), sigma = 0.5^toeplitz(1:4)))
X_Q = cbind(rep(1, n/10), rmvnorm(n/10, rep(0, 4), sigma = 0.5^toeplitz(1:4)))
beta = runif(5,-5,5)

y_P = c(X_P%*%beta + 5*X_P[,1]*X_P[,2] +rnorm(dim(X_P)[1]))
y_P.m = y_P
y_Q = c(X_Q%*%beta + 5*X_Q[,1]*X_Q[,2] +rnorm(dim(X_Q)[1]))
y_Q.m = y_Q



build_and_compile_model <- function(norm) {
  model <- keras_model_sequential() %>%
    norm() %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(1)
  
  model %>% compile(
    loss = 'mean_absolute_error',
    optimizer = optimizer_adam(0.001)
  )
  
  model
}

normalizer <- layer_normalization(axis = -1L)
normalizer %>% adapt(data = as.matrix(X_P))

dnn_model <- build_and_compile_model(normalizer)
history <- dnn_model %>% fit(
  as.matrix(X_P),
  as.matrix(y_P),
  validation_split = 0.2,
  verbose = 0,
  epochs = 1e2
)

test_results <- dnn_model %>% evaluate(
  as.matrix(X_Q),
  as.matrix(y_Q),
  verbose = 0
)

test_results 
TT2 =Sys.time()
mean((y_Q.m - X_Q%*%solve(t(X_P) %*% X_P) %*% t(X_P) %*% as.matrix(y_P.m))^2)
Sys.time() - TT1
Sys.time() - TT2

