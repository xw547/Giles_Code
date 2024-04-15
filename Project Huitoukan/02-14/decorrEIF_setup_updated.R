# Helper Functions

# Calculate the density ratio
density_ratio <- function(index, pw_index, pz_index, pwz) {
  curr_pw = pw_index$estimate[index]
  curr_pz = pz_index$estimate[index]
  curr_pwz = pwz$estimate[index]
  return(curr_pw * curr_pz / curr_pwz)
}

# density_ratio <- function(index, pw_index, pz_index, pwz) {
#   curr_pw = dnorm(pw_index$eval.points[index], 3, 1)
#   curr_pz = dnorm(pz_index$eval.points[index], 0.9, 1.09)
#   sigma_matrix <- matrix(c(1, 0.3, 0.3, 1.09), ncol = 2)
#   mean_vector <- c(3, 0.9)
#   curr_pwz = dmvnorm(pwz$eval.points[1,], mean_vector, sigma_matrix)
#   return(curr_pw * curr_pz / curr_pwz)
# }

# Calculate the first term 
first_term_z <- function(input_z, n = 100, pw, full_model, reduced_model) {
  samples = sample(pw$x, replace = T, size = n, prob = pw$y)
  predict_data_full = data.frame(w = samples, z = rep(input_z, n))
  predict_data_reduced = data.frame(z = rep(input_z, n))
  
  full_predictions = predict(full_model, predict_data_full)
  reduced_predictions = predict(reduced_model, predict_data_reduced)
  
  return(mean((full_predictions - reduced_predictions)^2))
}


# Calculate the second term 
second_term_w <- function(input_w, n = 100, pz, full_model, reduced_model) {
  samples = sample(pz$x, replace = T, size = n, prob = pz$y)
  predict_data_full = data.frame(w = rep(input_w, n), z = samples)
  predict_data_reduced = data.frame(z = samples)
  
  full_predictions = predict(full_model, predict_data_full)
  reduced_predictions = predict(reduced_model, predict_data_reduced)
  
  return(mean((full_predictions - reduced_predictions)^2))
}

# Calculate the third term 
third_term_wz <- function(index, input_w, input_z, input_y, pw_index, pz_index, pwz, full_model, reduced_model) {
  
  density = density_ratio(index, pw_index, pz_index, pwz)
  full_prediction = predict(full_model, data.frame(w = input_w, z = input_z))
  reduced_prediction = predict(reduced_model, data.frame(z = input_z))
  
  result = density * (full_prediction - reduced_prediction) * (input_y - full_prediction)
  
  return(result)
}


# Calculate psi_hat, the empirical estimator.
psi_hat <- function(n = 100, pw, pz, full_model, reduced_model) {
  sample_w <- sample(pw$x, size = n, replace = TRUE, prob = pw$y)
  sample_z <- sample(pz$x, size = n, replace = TRUE, prob = pz$y)
  predict_data_full <- data.frame(w = sample_w, z = sample_z)
  predict_data_reduced <- data.frame(z = sample_z)
  result <- mean((predict(full_model, predict_data_full) - predict(reduced_model, predict_data_reduced))^2)
  
  return(result)
}

DecorrLOCO <- function(index, inputw, inputz, inputy, pw, pw_index, 
                       pz, pz_index, pwz, full_model, reduced_model){
  third_term <- third_term_wz(index, inputw, inputz, inputy, pw, pw_index, pz, pz_index, pwz, full_model, reduced_model)
  first_term <- first_term_z(inputz, n = 100, pw, full_model, reduced_model)
  second_term <- second_term_w(inputw, n = 100, pz, full_model, reduced_model)
  return(2 * third_term + first_term + second_term)
}

EIF <- function(index, inputw, inputz, inputy, pw, pw_index, 
                pz, pz_index, pwz, full_model, reduced_model, psi_now){
  result = DecorrLOCO(index, inputw, inputz, inputy, pw, pw_index, 
                      pz, pz_index, pwz, full_model, reduced_model) - 2*psi_now
  return(result)
}


### Now, we proceed to the calculatin of \hat{\psi}
hat_psi <- function(w, z, y, pw, pw_index, 
                    pz, pz_index, pwz, full_model, reduced_model, psi_now){
  sample_size = length(y)
  phi_vec = c()
  phat_vec = c()
  for (index in 1:sample_size){
    phi_vec[index] <- EIF(index, w[index], z[index], y[index], pw, pw_index, 
                          pz, pz_index, pwz, full_model, reduced_model, psi_now)
    
  }
  phi_vec[abs(phi_vec)>100] = 0
  epsilon <- 0.05 * sign(c(mean(phi_vec)))
  second_term_hatpsi <- mean(exp(epsilon*phi_vec)*phi_vec)/mean(exp(epsilon*phi_vec))
  
  result <- mean(phi_vec) - second_term_hatpsi
  return(result)
}


# Calculate the likelihood (normalized update)
likelihood_nor <- function(gamma, y, pwz, tilde_phi) {
  likelihood_individual <- numeric(length(y))
  
  for (index in 1:length(y)) {
    likelihood_individual[index] <- pwz$estimate[index] * exp(gamma * tilde_phi[index]) - gamma * tilde_phi[index]
  }
  
  return(mean(likelihood_individual))
}

# Calculate the original likelihood
likelihood_ori <- function(gamma, y, pwz, eif_vec) {
  likelihood_individual <- numeric(length(y))
  
  for (index in 1:length(y)) {
    likelihood_individual[index] <- pwz$estimate[index] * exp(gamma * eif_vec[index]) - gamma * eif_vec[index]
  }
  
  return(mean(likelihood_individual))
}

likelihood_first <- function(epsilon, y, pycwz, eif_first){
  likelihood_individual <- numeric(length(y))
  
  for (index in 1:length(y)) {
    likelihood_individual[index] <- pycwz$condens[index] * exp(epsilon * eif_first[index]) - epsilon * eif_first[index]
  }
  
  return(mean(likelihood_individual))
}


