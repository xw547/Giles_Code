set.seed(2044)
n = 4000
n_1 = 2000
p = 10
rho = .5
beta = c(5, rep(0,9))
# beta = c(0, 5, rep(0,8))

sig = diag(rep(1,10))
sig[1,2] = rho
sig[2,1] = rho

X = rmvnorm(n, mean = rep(0,10), sigma = sig)
X_reduced = X[,-1]
Y = X%*%beta + rnorm(n, 0, 1)

full_data = data.frame(y = Y, X = X)
colnames(full_data) <- c("output", paste0("X", 1:10))
reduced_data = data.frame(y = Y, X_reduced = X_reduced)
colnames(reduced_data) <-  c("output", paste0("X", 2:10))

### Now we fit the models.
full_model = ranger(output~., data = full_data[1:n_1,])
reduced_model = ranger(output~., data= reduced_data[1:n_1,])


### We can now try to compare the naive estimator and the tMLE one-step one.
naive_LOCO = -mean((full_data[(n_1+1):n, 1]- predict(full_model, full_data[(n_1+1):n, -1])$predictions)^2 
                   - (full_data[(n_1+1):n, 1]- predict(reduced_model, reduced_data[(n_1+1):n, -1])$predictions)^2)


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
colnames(full_data_trunc) <- c("output", paste0("X", 1:10))
reduced_data_trunc = data.frame(y = Y_trunc, X_reduced = X_reduced[1:n_1,])
colnames(reduced_data_trunc) <- c("output", paste0("X", 2:10))

### Now we fit the models.
full_model_trunc = ranger(output~., data = full_data_trunc[1:n_1,])
reduced_model_trunc = ranger( output~., data= reduced_data_trunc[1:n_1,])


trunc_index = full_model_trunc$predictions>.005&full_model_trunc$predictions<(1-.005)
full_model_trunc$predictions = full_model_trunc$predictions[trunc_index]
Y_trunc = Y_trunc[trunc_index]


  eif_h = c()
  
  for (index in 1:length(Y_trunc)){
    eif_h[index] = first_term_noresidue_rf(index, full_data_trunc, 
                                           full_model_trunc, reduced_model_trunc)
  }
  
  regress_frame = data.frame(y = Y_trunc, curr_offset = qlogis(full_model_trunc$predictions), obs = eif_h)
  epsilon_model = glm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
  full_model_trunc$predictions = epsilon_model$fitted.values




epsilon_model$fitted.values = epsilon_model$fitted.values * dist + min



bb = mean((2*full_data[1:n_1, 1][trunc_index]
           - epsilon_model$fitted.values
           - reduced_model$predictions[trunc_index])*
            (epsilon_model$fitted.values - reduced_model$predictions[trunc_index]))

cc = mean((2*full_data_trunc[1:n_1, 1]
           - full_model_trunc$predictions
           - reduced_model_trunc$predictions)*
            (full_model_trunc$predictions - reduced_model_trunc$predictions))*dist^2

c(aa,bb,cc, naive_LOCO)
