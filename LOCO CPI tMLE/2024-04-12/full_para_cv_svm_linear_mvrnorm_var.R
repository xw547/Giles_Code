library(doParallel)
library(foreach)
library(mvtnorm)
library(gam)
library(MASS)
registerDoParallel(80)

source("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-12/helpers.R")
linear_truth = 25*(1-seq(.1,.9,.1)^2)

output = foreach(rho = seq(.1, .9, .1), .combine = "rbind")%do%{
  results = foreach(s = 1:200, .combine = "cbind")%dopar%{
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
    Y = X%*%beta + rnorm(n)
    
    
    
    full_data = data.frame(y = Y, X = X)
    colnames(full_data) <- c("output", 1:10)
    reduced_data = data.frame(y = Y, X_reduced = X_reduced)
    colnames(reduced_data) <- c("output", 2:10)
    
    
    X_2 = rmvnorm(n, mean = rep(0,10), sigma = sig)
    X_reduced_2 = X_2[,-1]
    Y_2 = X_2%*%beta + rnorm(n)
    full_data_2 = data.frame(y = Y_2, X = X_2)
    colnames(full_data_2) <- c("output", 1:10)
    reduced_data_2 = data.frame(y = Y_2, X_reduced_2 = X_reduced_2)
    colnames(reduced_data_2) <- c("output", 2:10)
    
    
    
    ### Now we fit the models.
    full_model = svm(output~., data = full_data[1:n_1,])
    reduced_model = svm(output~., data= reduced_data[1:n_1,])
    full_model_2 = svm(output~., data = full_data_2[1:n_1,])
    reduced_model_2 = svm(output~., data= reduced_data_2[1:n_1,])
    
    
    ### We can now try to compare the naive estimator and the tMLE one-step one.
    naive_LOCO = -((full_data_2[1:n_1, 1]- predict(full_model, full_data_2[1:n_1, -1]))^2 
                   - (full_data_2[1:n_1, 1]- predict(reduced_model, reduced_data_2[1:n_1, -1]))^2)
    naive_LOCO_2 = -((full_data[1:n_1, 1]- predict(full_model_2, full_data[1:n_1, -1]))^2 
                     - (full_data[1:n_1, 1]- predict(reduced_model_2, reduced_data[1:n_1, -1]))^2)
    
    
    ### Notice that the above two forms actually shares the same form, so the 
    ### One-step estimator isn't going to change anything.
    
    ### Now, we are ready to conduct the tMLE update.
    ### Notice that 
    ### Y|W, Z: 2(Y - \mu(W, Z))^2
    ### Y, W|Z: 0
    ### E : (\mu(W, Z) - \mu(Z))^2
    
    ### Data Processing
    Y_cv_1_1 <- rescaling(full_data[1:n_1,], 1)
    full_data_trunc_1_1 <- Y_cv_1_1[[1]]
    reduced_data_trunc_1_1 <- Y_cv_1_1[[2]]
    Y_trunc_1 <- Y_cv_1_1[[3]]
    dist_1 = Y_cv_1_1[[4]][1]
    min_1 = Y_cv_1_1[[4]][2]
    Y_cv_2_1 <- rescaling_cv(full_data_2[1:n_1,], 1, Y_cv_1_1[[4]][1], Y_cv_1_1[[4]][2])
    full_data_trunc_2_1 <- Y_cv_2_1[[1]]
    reduced_data_trunc_2_1 <- Y_cv_2_1[[2]]
    
    Y_cv_2_2 <- rescaling(full_data_2[1:n_1,], 1)
    full_data_trunc_2_2 <- Y_cv_2_2[[1]]
    reduced_data_trunc_2_2 <- Y_cv_2_2[[2]]
    Y_trunc_2 <- Y_cv_2_2[[3]]
    dist_2 = Y_cv_2_2[[4]][1]
    min_2  = Y_cv_2_2[[4]][2]
    Y_cv_1_2 <- rescaling_cv(full_data[1:n_1,], 1, Y_cv_2_2[[4]][1], Y_cv_2_2[[4]][2])
    full_data_trunc_1_2 <- Y_cv_1_2[[1]]
    reduced_data_trunc_1_2 <- Y_cv_1_2[[2]]
    
    full_model_trunc_1_1 <- svm(output~., data = full_data_trunc_1_1[1:n_1,])
    full_model_trunc_1_1_res <- svm(output~., data = full_data_trunc_1_1[1:n_1,])
    reduced_model_trunc_1_1 <-  svm(output~., data = reduced_data_trunc_1_1[1:n_1,])
    
    full_model_trunc_2_1_pred <- predict(full_model_trunc_1_1, full_data_trunc_2_1)
    reduced_model_trunc_2_1_pred <- predict(reduced_model_trunc_1_1, full_data_trunc_2_1)
    
    trunc_index_1_1 <- trunc(full_model_trunc_1_1$fitted.values, 
                             reduced_model_trunc_1_1$fitted.values)
    trunc_index_2_1 <- trunc(full_model_trunc_2_1_pred,  
                             reduced_model_trunc_2_1_pred)
    
    
    
    eif_h =  2*(full_model_trunc_1_1$fitted.values - reduced_model_trunc_1_1$fitted.values)[trunc_index_1_1]
    regress_frame = data.frame(y = Y_cv_1_1[[3]][trunc_index_1_1], 
                               curr_offset = qlogis(full_model_trunc_1_1$fitted.values[trunc_index_1_1]), 
                               obs = eif_h)
    
    
    epsilon_model = gsvm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
    full_model_trunc_1_1$fitted.values = epsilon_model$fitted.values
    
    
    
    full_model_trunc_2_1_pred = full_model_trunc_2_1_pred[trunc_index_2_1]
    Y_trunc_2 = Y_trunc_2[trunc_index_2_1]
    
    validataion_frame = data.frame(y = Y_trunc_2, 
                                   curr_offset = qlogis(full_model_trunc_2_1_pred),
                                   obs = 2*(full_model_trunc_2_1_pred - reduced_model_trunc_2_1_pred[trunc_index_2_1]))
    
    epsilon_model$fitted.values = predict(epsilon_model, newdata = validataion_frame, type = "response")
    
    epsilon_model$fitted.values = epsilon_model$fitted.values * dist_1 + min_1
    
    
    bb_1 =  ((2*full_data_2[1:n_1, 1][trunc_index_2_1]
              - epsilon_model$fitted.values
              - predict(reduced_model, reduced_data_2[1:n_1, ])[trunc_index_2_1])*
               (epsilon_model$fitted.values - predict(reduced_model, reduced_data_2[1:n_1, ])[trunc_index_2_1]))
    
    
    
    full_model_trunc_2_2 <- svm(output~., data = full_data_trunc_2_2[1:n_1,])
    full_model_trunc_2_2_res <- svm(output~., data = full_data_trunc_2_2[1:n_1,])
    reduced_model_trunc_2_2 <-  svm(output~., data = reduced_data_trunc_2_2[1:n_1,])
    
    full_model_trunc_1_2_pred <- predict(full_model_trunc_2_2, full_data_trunc_1_2)
    reduced_model_trunc_1_2_pred <- predict(reduced_model_trunc_2_2, full_data_trunc_1_2)
    
    trunc_index_2_2 <- trunc(full_model_trunc_2_2$fitted.values, 
                             reduced_model_trunc_2_2$fitted.values)
    trunc_index_1_2 <- trunc(full_model_trunc_1_2_pred,  
                             reduced_model_trunc_1_2_pred)
    
    
    
    eif_h =  2*(full_model_trunc_2_2$fitted.values - reduced_model_trunc_2_2$fitted.values)[trunc_index_2_2]
    regress_frame = data.frame(y = Y_cv_2_2[[3]][trunc_index_2_2], 
                               curr_offset = qlogis(full_model_trunc_2_2$fitted.values[trunc_index_2_2]), 
                               obs = eif_h)
    
    
    epsilon_model = gsvm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
    
    
    
    full_model_trunc_1_2_pred = full_model_trunc_1_2_pred[trunc_index_1_2]
    Y_trunc_1 = Y_trunc_1[trunc_index_1_2]
    
    validataion_frame = data.frame(y = Y_trunc_1, 
                                   curr_offset = qlogis(full_model_trunc_1_2_pred),
                                   obs = 2*(full_model_trunc_1_2_pred - reduced_model_trunc_1_2_pred[trunc_index_1_2]))
    
    epsilon_model$fitted.values = predict(epsilon_model, newdata = validataion_frame, type = "response")
    
    epsilon_model$fitted.values = epsilon_model$fitted.values * dist_2 + min_2
    
    bb_2 =  ((2*full_data[1:n_1, 1][trunc_index_1_2]
              - epsilon_model$fitted.values
              - predict(reduced_model_2, reduced_data[1:n_1, ])[trunc_index_1_2])*
               (epsilon_model$fitted.values - predict(reduced_model_2, reduced_data[1:n_1, ])[trunc_index_1_2]))
    
    
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
    full_model_trunc = svm(output~., data = full_data_trunc[1:n_1,])
    full_model_trunc_est = svm(output~., data = full_data_trunc[1:n_1,])
    reduced_model_trunc = svm( output~., data= reduced_data_trunc[1:n_1,])
    
    
    trunc_index = full_model_trunc$fitted.values>.005&full_model_trunc$fitted.values<(1-.005)
    full_model_trunc$fitted.values = full_model_trunc$fitted.values[trunc_index]
    Y_trunc = Y_trunc[trunc_index]
    
    
    eif_h = c(full_model_trunc$fitted.values - reduced_model_trunc$fitted.values[trunc_index])
    
    regress_frame = data.frame(y = Y_trunc, curr_offset = qlogis(full_model_trunc$fitted.values), obs = eif_h)
    epsilon_model = gsvm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
    
    
    
    epsilon_model$fitted.values = epsilon_model$fitted.values * dist + min
    
    aa_1 =  ((2*full_data[1:n_1, 1][trunc_index]
              - epsilon_model$fitted.values
              - reduced_model$fitted.values[trunc_index])*
               (epsilon_model$fitted.values - reduced_model$fitted.values[trunc_index]))
    dd_1 =  ((2*full_data_trunc[1:n_1, 1][trunc_index]
              - full_model_trunc_est$fitted.values[trunc_index]
              - reduced_model_trunc$fitted.values[trunc_index])*
               (full_model_trunc_est$fitted.values[trunc_index] - reduced_model_trunc$fitted.values[trunc_index]))*dist^2
    ### Bounded Logistic 
    y_2 = Y_2[1:n_1]
    dist_2 = (max(y_2) - min(y_2))
    min_2 = min(y_2)
    Y_trunc_2 = (y_2 - min(y_2))/(max(y_2) - min(y_2))
    full_data_trunc_2 = data.frame(y = Y_trunc_2, X = X_2[1:n_1,])
    colnames(full_data_trunc_2) <- c("output", paste0("X", 1:10))
    reduced_data_trunc_2 = data.frame(y = Y_trunc_2, X_reduced = X_reduced_2[1:n_1,])
    colnames(reduced_data_trunc_2) <- c("output", paste0("X", 2:10))
    
    ### Now we fit the models.
    full_model_trunc_2 = svm(output~., data = full_data_trunc_2[1:n_1,])
    full_model_trunc_est_2 = svm(output~., data = full_data_trunc_2[1:n_1,])
    reduced_model_trunc_2 = svm( output~., data= reduced_data_trunc_2[1:n_1,])
    
    
    
    
    trunc_index = full_model_trunc_2$fitted.values>.005&full_model_trunc_2$fitted.values<(1-.005)
    full_model_trunc_2$fitted.values = full_model_trunc_2$fitted.values[trunc_index]
    Y_trunc_2 = Y_trunc_2[trunc_index]
    
    
    eif_h = c(full_model_trunc_2$fitted.values - reduced_model_trunc_2$fitted.values[trunc_index])
    
    
    regress_frame = data.frame(y = Y_trunc_2, curr_offset = qlogis(full_model_trunc_2$fitted.values), obs = eif_h)
    epsilon_model = gsvm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")
    
    epsilon_model$fitted.values = epsilon_model$fitted.values * dist_2 + min_2
    
    
    
    
    aa_2 =  ((2*full_data_2[1:n_1, 1][trunc_index]
              - epsilon_model$fitted.values
              - reduced_model_2$fitted.values[trunc_index])*
               (epsilon_model$fitted.values - reduced_model_2$fitted.values[trunc_index]))
    
    dd_2 =   ((2*full_data_trunc_2[1:n_1, 1][trunc_index]
               - full_model_trunc_est_2$fitted.values[trunc_index]
               - reduced_model_trunc_2$fitted.values[trunc_index])*
                (full_model_trunc_est_2$fitted.values[trunc_index] - reduced_model_trunc_2$fitted.values[trunc_index]))*dist_2^2
    
    
    
    
    bb = mean(bb_1 + bb_2)/2
    bb_var = mean((bb_1 -linear_truth[i])^2 + (bb_2 - linear_truth[i])^2)/2
    
    dd = mean(dd_1 + dd_2)/2
    dd_var = mean((dd_1 -linear_truth[i])^2 + (dd_2 - linear_truth[i])^2)/2
    
    ee = mean(naive_LOCO + naive_LOCO_2)/2
    ee_var = mean((naive_LOCO -linear_truth[i])^2 + (naive_LOCO_2 - linear_truth[i])^2)/2
    
    c(bb, bb_var, dd, dd_var, ee, ee_var)}
  results
}

write.csv(output, file = "~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-12/full_svm_linear_var_1000.csv")


