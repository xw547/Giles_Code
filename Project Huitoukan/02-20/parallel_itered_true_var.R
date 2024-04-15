
  ### Data generation: In principle, we should be using the Friedman equation 
  ### to estimate the stuff. 
  s = 1
  set.seed(s+2013)
  
  ### Data generation: In principle, we should be using the Friedman equation 
  ### to estimate the stuff. 
  rho = 0.3
  
  w = rnorm(1000, 3, 1)
  y = 2 *w + rnorm(length(w), 0, 1)
  z = rho * w + rnorm(length(w), 0, 1)
  x = cbind(w,z)
  
  full_model = lm(y~1+w+z)
  reduced_model = lm(y~1+z)
  
  
  ### The following kde density estimates were used to estimate the density ratio 
  ### estimation. However, one thing to mention is that we chose to use index so 
  ### we may have a better speed. 
  
  pw_index = kde(w, eval.points = w, density = T)
  pz_index = kde(z, eval.points = z, density = T)
  pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)
  
  ### Cross-fitting to apply the normalization.
  ### Data generation: In principle, we should be using the Friedman equation 
  ### to estimate the stuff. 
  
  set.seed(s+2014)
  rho = 0.3
  
  w_2 = rnorm(1000, 3, 1)
  y_2 = 2 *w_2 + rnorm(length(w_2), 0, 1)
  z_2 = rho * w_2 + rnorm(length(w_2), 0, 1)
  x_2 = as.data.frame(cbind(w_2,z_2))
  colnames(x_2) <- c("w","z")
  
  full_model_2 = lm(y_2~1+w+z,data = x_2)
  reduced_model_2 = lm(y_2~1+z, data = x_2)
  
  
  ### The following kde density estimates were used to estimate the density ratio 
  ### estimation. However, one thing to mention is that we chose to use index so 
  ### we may have a better speed. 
  
  pw_index_2 = kde(w_2, eval.points = w_2, density = T)
  pw_index_2_1 = kde(w_2, eval.points = w, density = T)
  pz_index_2 = kde(z_2, eval.points = z_2, density = T)
  pz_index_2_1 = kde(z_2, eval.points = z, density = T)
  pwz_2 = kde(cbind(w_2,z_2), eval.points = cbind(w_2,z_2), density = T)
  pwz_2_1 = kde(cbind(w_2,z_2), eval.points = cbind(w,z), density = T)
  
  
  decorr_eif = c()
  eif_inital = c()
  
  
  tol = 1e-5
  counter = 1 
  eps = 1
  psi_now = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)
  
  if(eps == 1){
    for (index in seq_along(y)) {
      decorr_eif[index] = DecorrLOCO(index, w[index], z[index], y[index], pw_index, pw_index, pz_index, pz_index, pwz, full_model, reduced_model)
      eif_inital[index] = EIF(index, w[index], z[index], y[index], pw_index, pw_index, pz_index, pz_index, pwz, full_model, reduced_model, psi_now)
    }
    eif_estimate = 1/2*mean(decorr_eif)
    var_initial_est = mean(eif_inital^2)}
  
  while(eps > tol && counter <50){
    
    psi_now = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)
    
    eif_vec = c()
    
    for (index in 1:length(y)){
      eif_vec[index] = EIF(index, w[index], z[index], y[index], pw_index, pw_index, 
                           pz_index, pz_index, pwz, full_model, reduced_model, psi_now)
    }
    
    epsilon = optimise(likelihood_ori, c(-1/2, 1/2), pwz = pwz, eif_vec = eif_vec, y = y)
    
    ### Updating the distributions
    hat_phi_applied = exp(eif_vec * epsilon$objective)
    pw_index$estimate = (pw_index$estimate * hat_phi_applied) / sum(pw_index$estimate * hat_phi_applied)
    pz_index$estimate = (pz_index$estimate * hat_phi_applied) / sum(pz_index$estimate * hat_phi_applied)
    pwz$estimate = (pwz$estimate * hat_phi_applied) / sum(pwz$estimate * hat_phi_applied)
    counter = counter+1
    eps = epsilon$objective
  }
  
  
  eif_vec = c()
  decorr_vec = c()
  
  
  for (index in seq_along(y)) {
    psi_now = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)
    eif_vec[index] = EIF(index, w[index], z[index], y[index], pw_index, pw_index, pz_index, pz_index, pwz, full_model, reduced_model, psi_now)
    decorr_vec[index] = DecorrLOCO(index, w[index], z[index], y[index], pw_index, pw_index, pz_index, pz_index, pwz, full_model, reduced_model)
  }
  
  c(eif_estimate, 1/2 * mean(decorr_vec, na.rm = T), var_initial_est, mean(eif_vec^2, na.rm = T))
  
  

  set.seed(s+2013)
  
  ### Data generation: In principle, we should be using the Friedman equation 
  ### to estimate the stuff. 
  rho = 0.3
  
  w = rnorm(1000, 3, 1)
  y = 2 *w + rnorm(length(w), 0, 1)
  z = rho * w + rnorm(length(w), 0, 1)
  x = cbind(w,z)
  
  full_model = lm(y~1+w+z)
  reduced_model = lm(y~1+z)
  
  
  ### The following kde density estimates were used to estimate the density ratio 
  ### estimation. However, one thing to mention is that we chose to use index so 
  ### we may have a better speed. 
  
  pw_index = kde(w, eval.points = w, density = T)
  pz_index = kde(z, eval.points = z, density = T)
  pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)
  
  ### Cross-fitting to apply the normalization.
  ### Data generation: In principle, we should be using the Friedman equation 
  ### to estimate the stuff. 
  
  set.seed(s+2014)
  rho = 0.3
  
  w_2 = rnorm(1000, 3, 1)
  y_2 = 2 *w_2 + rnorm(length(w_2), 0, 1)
  z_2 = rho * w_2 + rnorm(length(w_2), 0, 1)
  x_2 = as.data.frame(cbind(w_2,z_2))
  colnames(x_2) <- c("w","z")
  
  full_model_2 = lm(y_2~1+w+z,data = x_2)
  reduced_model_2 = lm(y_2~1+z, data = x_2)
  
  
  ### The following kde density estimates were used to estimate the density ratio 
  ### estimation. However, one thing to mention is that we chose to use index so 
  ### we may have a better speed. 
  
  pw_index_2 = kde(w_2, eval.points = w_2, density = T)
  pw_index_2_1 = kde(w_2, eval.points = w, density = T)
  pz_index_2 = kde(z_2, eval.points = z_2, density = T)
  pz_index_2_1 = kde(z_2, eval.points = z, density = T)
  pwz_2 = kde(cbind(w_2,z_2), eval.points = cbind(w_2,z_2), density = T)
  pwz_2_1 = kde(cbind(w_2,z_2), eval.points = cbind(w,z), density = T)
  
  psi_final = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)
  
  ### gamma
  tol = 1e-5
  counter = 1 
  gma = 1
  
  while(gma > tol && counter <2){
    psi_now = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)
    
    psi_curr_1 <- hat_psi(w, z, y, pw_index, pw_index, 
                          pz_index, pz_index, pwz, full_model, reduced_model, psi_now)
    psi_curr_2_1 <- hat_psi(w, z, y, pw_index_2, pw_index_2_1, 
                            pz_index_2, pz_index_2_1, pwz_2_1, full_model_2, reduced_model_2, psi_now)
    phi_curr = c(psi_curr_1, psi_curr_2_1 )
    
    tilde_phi <- c()
    eif_vec = c()
    
    
    for (index in 1:length(y)){
      tilde_phi[index] = 
        (EIF(index, w[index], z[index], y[index], pw_index, pw_index, 
             pz_index, pz_index, pwz, full_model, reduced_model, psi_now) * phi_curr[1] +
           EIF(index, w[index], z[index], y[index], pw_index_2_1, pw_index_2_1, 
               pz_index_2_1, pz_index_2_1, pwz_2_1, full_model_2, reduced_model_2, psi_now) * phi_curr[2])/
        norm(phi_curr, "2")
    }
    
    gamma = optimise(likelihood_nor, c(-1/2, 1/2), pwz = pwz, tilde_phi = tilde_phi, y = y)
    
    tilde_phi_applied = exp(tilde_phi * gamma$objective)
    pw_index$estimate = (pw_index$estimate * tilde_phi_applied) / sum(pw_index$estimate * tilde_phi_applied)
    pz_index$estimate = (pz_index$estimate * tilde_phi_applied) / sum(pz_index$estimate * tilde_phi_applied)
    pwz$estimate = (pwz$estimate * tilde_phi_applied) / sum(pwz$estimate * tilde_phi_applied)
    gma = gamma$objective
    counter = counter +1
  }
  
  eif_vec_ori = c()
  decorr_vec_ori = c()
  
  for (index in seq_along(y)) {
    psi_now = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)
    eif_vec_ori[index] = EIF(index, w[index], z[index], y[index], pw_index, pw_index, pz_index, pz_index, pwz, full_model, reduced_model, psi_now)
    decorr_vec_ori[index] = DecorrLOCO(index, w[index], z[index], y[index], pw_index, pw_index, pz_index, pz_index, pwz, full_model, reduced_model)
  }
  
  
  # sample_w = rnorm(1000, 3, 1)
  # sample_z = rnorm(1000, 0.9, 1.09)
  # 
  # predict_data_full = data.frame(x = cbind(sample_w, sample_z))
  # colnames(predict_data_full) = c("w", "z")
  # predict_data_reduced = data.frame(z = sample_z)
  
  
  psi_now = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)
  c(1/2 * mean(decorr_vec_ori, na.rm = T), psi_final, mean(eif_vec_ori^2, na.rm = T)) 
  #mean((predict(full_model, predict_data_full)-
  #predict(reduced_model, predict_data_reduced))^2)
