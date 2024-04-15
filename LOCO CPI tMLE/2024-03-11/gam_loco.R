### In this code I'll try to implement the tmle for LOCO estimates.
### This shall include a two cases, the one step estimator and the 
### tMLE version. 

### Our setting shall largely follows the paper from Williamson JASA since we 
### Shares a similar setting.

### Libraries used

library(MASS)
library(gam)

### Data Generating Process
set.seed(2025)
n = 1500
p = 10

beta = c((5), rep(0,9))
X = matrix(rnorm(n*p), ncol = 10)
X_reduced = X[,-1]
Y = c(X%*%beta) + rnorm(n, 1, 2)

full_data = data.frame(y = Y, X = X)
colnames(full_data) <- c("output", 1:10)
reduced_data = data.frame(y = Y, X_reduced = X_reduced)
colnames(reduced_data) <- c("output", 2:10)

### Now we fit the models.
full_model = gam(output~., data = full_data[1:1000,],  family = "gaussian")
reduced_model = gam(data= reduced_data[1:1000,], output~., family = "gaussian")


### We can now try to compare the naive estimator and the tMLE one-step one.
naive_LOCO = -mean((full_data[1001:1500, 1]- predict(full_model, full_data[1001:1500, -1]))^2 
                   - (full_data[1001:1500, 1]- predict(reduced_model, reduced_data[1001:1500, -1]))^2)
OneStep_LOCO = mean((2*full_data[1001:1500, 1] 
                     - predict(full_model, full_data[1001:1500, -1]) 
                     - predict(reduced_model, reduced_data[1001:1500, -1]))*
                      (predict(full_model, full_data[1001:1500, -1]) -
                         predict(reduced_model, reduced_data[1001:1500, -1])))


### Notice that the above two forms actually shares the same form, so the 
### One-step estimator isn't going to change anything.

### Now, we are ready to conduct the tMLE update.
### Notice that 
### Y|W, Z: 2(Y - \mu(W, Z))^2
### Y, W|Z: 0
### E : (\mu(W, Z) - \mu(Z))^2


### Bounded Logistic
y = Y[1:1000]
dist = (max(y) - min(y))
min = min(y)
Y_trunc = (y - min(y))/(max(y) - min(y))
full_data_trunc = data.frame(y = Y_trunc, X = X[1:1000,])
colnames(full_data_trunc) <- c("output", 1:10)
reduced_data_trunc = data.frame(y = Y_trunc, X_reduced = X_reduced[1:1000,])
colnames(reduced_data_trunc) <- c("output", 2:10)

### Now we fit the models.
full_model_trunc = glm(output~., data = full_data_trunc[1:1000,],  family = "gaussian")
reduced_model_trunc = glm( output~., data= reduced_data_trunc[1:1000,], family = "gaussian")


trunc_index = full_model_trunc$fitted.values>.005&full_model_trunc$fitted.values<(1-.005)
full_model_trunc$fitted.values = full_model_trunc$fitted.values[trunc_index]
Y_trunc = Y_trunc[trunc_index]


eif_h = c()

for (index in 1:length(Y_trunc)){
  eif_h[index] = first_term_noresidue(index, full_data_trunc, full_model_trunc)
}

regress_frame = data.frame(y = Y_trunc, curr_offset = qlogis(full_model_trunc$fitted.values), obs = eif_h)
epsilon_model = glm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
full_model_trunc$fitted.values = epsilon_model$fitted.values



epsilon = coef(epsilon_model)
coef(epsilon_model)

epsilon_model$fitted.values = epsilon_model$fitted.values * dist + min



mean((2*full_data[1:1000, 1][trunc_index]
      - epsilon_model$fitted.values
      - reduced_model$fitted.values[trunc_index])*
       (epsilon_model$fitted.values - reduced_model$fitted.values[trunc_index]))

mean((2*full_data[1:1000, 1]
      - full_model$fitted.values
      - reduced_model$fitted.values)*
       (full_model$fitted.values - reduced_model$fitted.values))


### Linear Update

# 
#   eif_h = c()
# 
#   for (index in 1:1000){
#     eif_h[index] = first_term(index, full_data, full_model)
#   }
# 
#   regress_frame = data.frame(y = Y[1:1000], curr_offset = (full_model$fitted.values), obs = eif_h)
#   epsilon_model <- glm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame, family = "gaussian")
# 
#   full_model$fitted.values = epsilon_model$fitted.values
#   epsilon = abs(epsilon_model$coefficients)
#   epsilon



### Results

mean((2*full_data[1:1000, 1]
      - epsilon_model$fitted.values
      - reduced_model$fitted.values)*
       (epsilon_model$fitted.values - reduced_model$fitted.values))

# mean((2*full_data[1:1000, 1]
#       - full_model$fitted.values
#       - reduced_model$fitted.values)*
#        (full_model$fitted.values - reduced_model$fitted.values))



### Helper Functions
first_term_noresidue <-
  function(index, full_data, full_model) {
    mean(2*(full_data[index, 1] 
            - full_model$fitted.values[index]))
  }

first_term <-
  function(index, full_data, full_model) {
    mean(2*(full_data[index, 1] 
            - full_model$fitted.values[index])^2)
  }


logit <- function(p) {
  log(p / (1 - p))
}

