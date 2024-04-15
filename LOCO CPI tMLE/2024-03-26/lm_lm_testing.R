set.seed(s + 20044)
n = 1500
n_1 = 1000
p = 10

beta = c(5, rep(0,9))

sig = diag(rep(1,10))
sig[1,2] = rho
sig[2,1] = rho

X = rmvnorm(n, mean = rep(0,10), sigma = sig)
X_reduced = X[,-1]
Y = X%*%beta + rnorm(n, 0, 1)

full_data = data.frame(y = Y, X = X)
colnames(full_data) <- c("output", 1:10)
reduced_data = data.frame(y = Y, X_reduced = X_reduced)
colnames(reduced_data) <- c("output", 2:10)

### Now we fit the models.
full_model = gam(output~., data = full_data[1:n_1,],  family = "gaussian")
reduced_model = gam(data= reduced_data[1:n_1,], output~., family = "gaussian")


### We can now try to compare the naive estimator and the tMLE one-step one.
naive_LOCO = -mean((full_data[(n_1+1):n, 1]- predict(full_model, full_data[(n_1+1):n, -1]))^2 
                   - (full_data[(n_1+1):n, 1]- predict(reduced_model, reduced_data[(n_1+1):n, -1]))^2)


### Notice that the above two forms actually shares the same form, so the 
### One-step estimator isn't going to change anything.

### Now, we are ready to conduct the tMLE update.
### Notice that 
### Y|W, Z: 2(Y - \mu(W, Z))^2
### Y, W|Z: 0
### E : (\mu(W, Z) - \mu(Z))^2


### Bounded Logistic
y = Y[1:n_1]
dist = (max(y) - min(y))
min = min(y)
Y_trunc = (y - min(y))/(max(y) - min(y))
full_data_trunc = data.frame(y = Y_trunc, X = X[1:n_1,])
colnames(full_data_trunc) <- c("output", 1:10)
reduced_data_trunc = data.frame(y = Y_trunc, X_reduced = X_reduced[1:n_1,])
colnames(reduced_data_trunc) <- c("output", 2:10)

### Now we fit the models.
full_model_trunc = gam(output~., data = full_data_trunc[1:n_1,],  family = "gaussian")
full_model_trunc_est = gam(output~., data = full_data_trunc[1:n_1,],  family = "gaussian")
reduced_model_trunc = gam( output~., data= reduced_data_trunc[1:n_1,], family = "gaussian")


trunc_index = full_model_trunc$fitted.values>.005&full_model_trunc$fitted.values<(1-.005)
full_model_trunc$fitted.values = full_model_trunc$fitted.values[trunc_index]
Y_trunc = Y_trunc[trunc_index]

epsilon = 1
iter = 1
while(epsilon > 1e-4&iter<200){
  eif_h = c(full_model_trunc$fitted.values - reduced_model_trunc$fitted.values[trunc_index])
  
  
  regress_frame = data.frame(y = Y_trunc, curr_offset = qlogis(full_model_trunc$fitted.values), obs = eif_h)
  epsilon_model = glm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
  full_model_trunc$fitted.values = epsilon_model$fitted.values
  
  
  
  epsilon = coef(epsilon_model)
  iter = iter+1
}
aa = coef(epsilon_model)

epsilon_model$fitted.values = epsilon_model$fitted.values * dist + min



bb = mean((2*full_data[1:n_1, 1][trunc_index]
           - epsilon_model$fitted.values
           - reduced_model$fitted.values[trunc_index])*
            (epsilon_model$fitted.values - reduced_model$fitted.values[trunc_index]))

cc = mean((2*full_data[1:n_1, 1]
           - full_model$fitted.values
           - reduced_model$fitted.values)*
            (full_model$fitted.values - reduced_model$fitted.values))

dd = mean((2*full_data_trunc[1:n_1, 1][trunc_index]
           - full_model_trunc_est$fitted.values[trunc_index]
           - reduced_model_trunc$fitted.values[trunc_index])*
            (full_model_trunc_est$fitted.values[trunc_index] - reduced_model_trunc$fitted.values[trunc_index]))*dist^2


c(bb,cc, dd, naive_LOCO)