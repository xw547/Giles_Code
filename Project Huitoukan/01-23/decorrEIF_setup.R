# Helper Functions

# Calculate the density ratio
density_ratio <- function(index, pw_index, pz_index, pwz) {
  curr_pw = pw_index$estimate[index]
  curr_pz = pz_index$estimate[index]
  curr_pwz = pwz$estimate[index]
  return(curr_pw * curr_pz / curr_pwz)
}

# Calculate the first term 
first_term_z <- function(input_z, n = 100, pw, full_model, reduced_model) {
  samples = sample(pw$x, size = n, prob = pw$y)
  predict_data_full = data.frame(w = samples, z = rep(input_z, n))
  predict_data_reduced = data.frame(z = rep(input_z, n))
  
  full_predictions = predict(full_model, predict_data_full)
  reduced_predictions = predict(reduced_model, predict_data_reduced)
  
  return(mean((full_predictions - reduced_predictions)^2))
}

# Calculate the second term 
second_term_w <- function(input_w, n = 100, pz, full_model, reduced_model) {
  samples = sample(pz$x, size = n, prob = pz$y)
  predict_data_full = data.frame(w = rep(input_w, n), z = samples)
  predict_data_reduced = data.frame(z = samples)
  
  full_predictions = predict(full_model, predict_data_full)
  reduced_predictions = predict(reduced_model, predict_data_reduced)
  
  return(mean((full_predictions - reduced_predictions)^2))
}

# Calculate the third term 
third_term_wz <- function(index, input_w, input_z, input_y, pw, pw_index, pz, pz_index, pwz, full_model, reduced_model) {
  density = density_ratio(index, pw_index, pz_index, pwz)
  full_prediction = predict(full_model, data.frame(w = input_w, z = input_z))
  reduced_prediction = predict(reduced_model, data.frame(w = input_w))
  
  result = density * (full_prediction - reduced_prediction) * (input_y - full_prediction)
  
  return(result)
}


psi_hat <- function(n = 100, pw, pz, full_model, reduced_model){
  sample_w = sample(pw$x, size = n, prob = pw$y)
  sample_z = sample(pz$x, size = n, prob = pz$y)
  predict_data_full = data.frame(x = cbind(sample_w, sample_z))
  colnames(predict_data_full) = c("w", "z")
  predict_data_reduced = data.frame(z = sample_z)
  result = mean((predict(full_model, predict_data_full)-
                   predict(reduced_model, predict_data_reduced))^2)
  return(result)
}

DecorrLOCO <- function(index, inputw, inputz, inputy, pw, pw_index, 
                       pz, pz_index, pwz, full_model, reduced_model){
  third_term = third_term_wz(index, inputw, inputz, inputy, pw, pw_index, 
                             pz, pz_index, pwz, full_model, reduced_model)
  first_term = first_term_z(inputz, n = 100, pw, full_model, reduced_model)
  second_term = second_term_w(inputw, n = 100, pw, full_model, reduced_model)
  return(2*third_term+ first_term+ second_term)
}

EIF <- function(index, inputw, inputz, inputy, pw, pw_index, 
                pz, pz_index, pwz, full_model, reduced_model, psi_now){
  result = DecorrLOCO(index, inputw, inputz, inputy, pw, pw_index, 
                      pz, pz_index, pwz, full_model, reduced_model) - 2*psi_now
  return(result)
}


### Now, we proceed to the calculatin of \hat{\psi}
hat_psi <- function(w, z, y, pw, pw_index, 
                    pz, pz_index, pwz, full_model, reduced_model){
  sample_size = length(y)
  phi_vec = c()
  phat_vec = c()
  for (index in 1:sample_size){
    phi_vec[index] <- EIF(index, w[index], z[index], y[index], pw, pw_index, 
                          pz, pz_index, pwz, full_model, reduced_model, psi_hat(n = 100, pw, pz, full_model, reduced_model))
    phat_vec[index] <- pwz$estimate[index]
  }
  phi_vec[abs(phi_vec)>100] = 0
  #print(max(phi_vec))
  #print(min(phi_vec))
  #print(mean(phat_vec*exp(epsilon*phi_vec)*phi_vec))
  #print(mean(phat_vec*exp(epsilon*phi_vec)))
  
  result <-mean(phi_vec) - mean(phat_vec*phi_vec)
  return(result)
}


likelihood_nor <- function(gamma){
  likelihood_individual = c()
  for (index in 1:length(y)){
    likelihood_individual[index] <- pwz$estimate[index]*exp(gamma*tilde_phi[index]) 
    - gamma*tilde_phi[index]
  }
  return(mean(likelihood_individual))
}


### We can start with the revised likelihood function and try to optimize over
### it. 

likelihood_ori <- function(gamma){
  likelihood_individual = c()
  for (index in 1:length(y)){
    likelihood_individual[index] <- pwz$estimate[index]*exp(gamma*eif_vec[index]) 
    - gamma*eif_vec[index]
  }
  return(mean(likelihood_individual))
}