library(torch)
library(mvtnorm)
# Create a synthetic dataset
set.seed(42) # For reproducibility
n = 5*1e3

X_P = cbind(rep(1, n), rmvnorm(n, rep(0, 4), sigma = 0.5^toeplitz(1:4)))
X_Q = cbind(rep(1, n/10), rmvnorm(n/10, rep(0, 4), sigma = 0.5^toeplitz(1:4)))
beta = c(1:5)

y_P = c(X_P%*%beta +rnorm(dim(X_P)[1]))
y_P.m = y_P
y_Q = c(X_Q%*%beta +rnorm(dim(X_Q)[1]))
y_Q.m = y_Q

y_P = torch_tensor(y_P, dtype = torch_float32())$view(c(-1, 1))
y_Q = torch_tensor(y_Q, dtype = torch_float32())$view(c(-1, 1))

model <- nn_module(
  initialize = function() {
    self$hidden <- nn_linear(in_features = n_features, out_features = 50)
    self$output <- nn_linear(in_features = 50, out_features = 1)
  },
  forward = function(x) {
    x %>% 
      self$hidden() %>%
      nnf_relu() %>%
      self$output()
  }
)

n_features <- ncol(X_P) # Assuming x_train is a matrix or data frame
model <- model()
loss_fn <- nnf_mse_loss
optimizer <- optim_adam(model$parameters)


for (epoch in 1:1e3) { # You can adjust the number of epochs
  optimizer$zero_grad()
  output <- model(X_P)
  loss <- loss_fn(output, y_P)
  loss$backward()
  optimizer$step()
}


predictions <- aa(X_Q)
test_loss <- mean(loss_fn(predictions, y_Q)$item())
cat(sprintf("Test loss: %f\n", test_loss))


mean((y_Q.m - X_Q%*%solve(t(X_P) %*% X_P) %*% t(X_P) %*% as.matrix(y_P.m))^2)




