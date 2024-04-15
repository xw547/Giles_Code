library(doParallel)
library(foreach)
library(ranger)
library(mvtnorm)
library(mgcv)
library(MASS)
registerDoParallel(80)

source("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-05/helpers.R")




rho = .5
set.seed(20044)
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
colnames(full_data) <- c("output", paste0("X", 1:10))
reduced_data = data.frame(y = Y, X_reduced = X_reduced)
colnames(reduced_data) <- c("output", paste0("X", 2:10))


X_2 = as.matrix(cbind(expand.grid(x1 = seq(-2, 2, length.out = 100), x2 = seq(-2, 2, length.out = 100),
                                  x3 = seq(-2, 2, length.out = 100)), matrix(rep(0, (p-3)*1e6), ncol = (p-3)))) 
X_reduced_2 = X_2[,-1]
Y_2 = X_2%*%beta + rnorm(1e6, 0, 1)
full_data_2 = data.frame(y = Y_2, X = X_2)
colnames(full_data_2) <- c("output", paste0("X", 1:10))
reduced_data_2 = data.frame(y = Y_2, X_reduced_2 = X_reduced_2)
colnames(reduced_data_2) <- c("output",  paste0("X", 2:10))



### Now we fit the models.
full_rf = ranger(output~., data = full_data[1:n_1,])
reduced_rf = ranger(output~., data= reduced_data[1:n_1,])
full_rf_2 = ranger(output~., data = full_data_2[1:n_1,])
reduced_rf_2 = ranger(output~., data= reduced_data_2[1:n_1,])



Y_cv_1_1 <- rescaling(full_data[1:n_1,], 1)
full_data_trunc_1_1 <- Y_cv_1_1[[1]]
reduced_data_trunc_1_1 <- Y_cv_1_1[[2]]
Y_trunc_1 <- Y_cv_1_1[[3]]
dist_1 = Y_cv_1_1[[4]][1]
min_1 = Y_cv_1_1[[4]][2]
Y_cv_2_1 <- rescaling_cv(full_data_2[,], 1, Y_cv_1_1[[4]][1], Y_cv_1_1[[4]][2])
full_data_trunc_2_1 <- Y_cv_2_1[[1]]
reduced_data_trunc_2_1 <- Y_cv_2_1[[2]]

Y_cv_2_2 <- rescaling(full_data_2[1:n_1,], 1)
full_data_trunc_2_2 <- Y_cv_2_2[[1]]
reduced_data_trunc_2_2 <- Y_cv_2_2[[2]]
Y_trunc_2 <- Y_cv_2_2[[3]]



full_rf_trunc_1_1 <- ranger(output~., data = full_data_trunc_1_1[1:n_1,])
full_rf_trunc_1_1_res <- ranger(output~., data = full_data_trunc_1_1[1:n_1,])
reduced_rf_trunc_1_1 <-  ranger(output~., data = reduced_data_trunc_1_1[1:n_1,])

full_rf_trunc_2_1_pred <- predict(full_rf_trunc_1_1, full_data_trunc_2_1)$predictions
reduced_rf_trunc_2_1_pred <- predict(reduced_rf_trunc_1_1, full_data_trunc_2_1)$predictions

trunc_index_1_1 <- trunc(full_rf_trunc_1_1$predictions, 
                         reduced_rf_trunc_1_1$predictions)
trunc_index_2_1 <- trunc(full_rf_trunc_2_1_pred,  
                         reduced_rf_trunc_2_1_pred)



eif_h =  2*(full_rf_trunc_1_1$predictions - reduced_rf_trunc_1_1$predictions)[trunc_index_1_1]
regress_frame = data.frame(y = Y_cv_1_1[[3]][trunc_index_1_1], 
                           rf_offset = qlogis(full_rf_trunc_1_1$predictions[trunc_index_1_1]), 
                           obs = eif_h)


epsilon_rf = glm(y ~ 0 + offset(rf_offset) + obs, data = regress_frame,  family = "quasibinomial")
full_rf_trunc_1_1$predictions = epsilon_rf$fitted.values



full_rf_trunc_2_1_pred = full_rf_trunc_2_1_pred[trunc_index_2_1]
Y_trunc_2 = Y_trunc_2[trunc_index_2_1]


validataion_frame = data.frame(y = Y_trunc_2, 
                               rf_offset = qlogis(full_rf_trunc_2_1_pred),
                               obs = 2*(full_rf_trunc_2_1_pred - reduced_rf_trunc_2_1_pred[trunc_index_2_1]))

epsilon_rf$fitted.values = predict(epsilon_rf, newdata = validataion_frame, type = "response")

epsilon_rf$fitted.values = epsilon_rf$fitted.values * dist_1 + min_1

rf_reduced = predict(reduced_rf, reduced_data_2[, ])$predictions[trunc_index_2_1]

rf_1 = (2*full_data_2[, 1][trunc_index_2_1]
             - epsilon_rf$fitted.values
             - rf_reduced )*
              (epsilon_rf$fitted.values -rf_reduced)

################################################################
################################################################
################################################################
################################################################
################################################################


### Now we fit the models.
full_lm = lm(output~., data = full_data[1:n_1,])
reduced_lm = lm(output~., data= reduced_data[1:n_1,])




full_lm_trunc_1_1 <- lm(output~., data = full_data_trunc_1_1[1:n_1,])
full_lm_trunc_1_1_res <- lm(output~., data = full_data_trunc_1_1[1:n_1,])
reduced_lm_trunc_1_1 <-  lm(output~., data = reduced_data_trunc_1_1[1:n_1,])

full_lm_trunc_2_1_pred <- predict(full_lm_trunc_1_1, full_data_trunc_2_1)
reduced_lm_trunc_2_1_pred <- predict(reduced_lm_trunc_1_1, full_data_trunc_2_1)

trunc_index_1_1 <- trunc(full_lm_trunc_1_1$fitted.values, 
                         reduced_lm_trunc_1_1$fitted.values)
trunc_index_2_1 <- trunc(full_lm_trunc_2_1_pred,  
                         reduced_lm_trunc_2_1_pred)


eif_h =  2*(full_lm_trunc_1_1$fitted.values - reduced_lm_trunc_1_1$fitted.values)[trunc_index_1_1]
regress_frame = data.frame(y = Y_cv_1_1[[3]][trunc_index_1_1], 
                           curr_offset = qlogis(full_lm_trunc_1_1$fitted.values[trunc_index_1_1]), 
                           obs = eif_h)


epsilon_lm = glm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
full_lm_trunc_1_1$fitted.values = epsilon_lm$fitted.values



full_lm_trunc_2_1_pred = full_lm_trunc_2_1_pred[trunc_index_2_1]
Y_trunc_2 = Y_trunc_2[trunc_index_2_1]

validataion_frame = data.frame(y = Y_trunc_2, 
                               curr_offset = qlogis(full_lm_trunc_2_1_pred),
                               obs = 2*(full_lm_trunc_2_1_pred - reduced_lm_trunc_2_1_pred[trunc_index_2_1]))

epsilon_lm$fitted.values = predict(epsilon_lm, newdata = validataion_frame, type = "response")

epsilon_lm$fitted.values = epsilon_lm$fitted.values * dist_1 + min_1

lm_1 = (2*full_data_2[, 1][trunc_index_2_1]
             - epsilon_lm$fitted.values
             - predict(reduced_lm, reduced_data_2[, ])[trunc_index_2_1])*
              (epsilon_lm$fitted.values - predict(reduced_lm, reduced_data_2[, ])[trunc_index_2_1])

truth = (2*full_data_2[, 1][trunc_index_2_1]
         - 5*full_data_2[,2][trunc_index_2_1]
         - 5*rho*full_data_2[,3][trunc_index_2_1])*
           (5*full_data_2[,2][trunc_index_2_1] - 5*rho*full_data_2[,3][trunc_index_2_1])

mise_rf = mean((rf_1 - truth)^2)
mise_lm = mean((lm_1 - truth)^2)

c(mise_rf, mise_lm)
