library(ks)

set.seed(2013)
rho = 0.3

w = rnorm(1000, 3, 1)
y = 2 * w + rnorm(length(w), 0, 1)
z = rho * w + rnorm(length(w), 0, 1)
x = cbind(w,z)

full_model = lm(y ~ 1 + w + z)
reduced_model = lm(y ~ 1 + z)

### Density estimation:
### Notice that we are only interested in the density estimation of existing
### points.

pw = kde(w, density = T)
pz = kde(z, density = T)

### The following kde density estimates were used to estimate the density ratio 
### estimation. However, one thing to mention is that we chose to use index so 
### we may have a better speed. 

pw_index = kde(w, eval.points = w, density = T)
pz_index = kde(z, eval.points = z, density = T)
pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)

### Cross-fitting to apply the normalization.
### Data generation: In principle, we should be using the Friedman equation 
### to estimate the stuff. 

set.seed(2014)
rho = 0.3

w_2 = rnorm(1000, 3, 1)
y_2 = 2 *w_2 + rnorm(length(w_2), 0, 1)
z_2 = rho * w_2 + rnorm(length(w_2), 0, 1)
x_2 = as.data.frame(cbind(w_2,z_2))
colnames(x_2) <- c("w","z")

full_model_2 = lm(y_2~1+w+z,data = x_2)
reduced_model_2 = lm(y_2~1+z, data = x_2)

### Density estimation:
### Notice that we are only interested in the density estimation of existing
### points.

pw_2 = kde(w_2, density = T)
pz_2 = kde(z_2, density = T)

### The following kde density estimates were used to estimate the density ratio 
### estimation. However, one thing to mention is that we chose to use index so 
### we may have a better speed. 

pw_index_2 = kde(w_2, eval.points = w_2, density = T)
pw_index_2_1 = kde(w_2, eval.points = w, density = T)
pz_index_2 = kde(z_2, eval.points = z_2, density = T)
pz_index_2_1 = kde(z_2, eval.points = z, density = T)
pwz_2 = kde(cbind(w_2,z_2), eval.points = cbind(w_2,z_2), density = T)
pwz_2_1 = kde(cbind(w_2,z_2), eval.points = cbind(w,z), density = T)


### Simple simulation to see if our result works
index = 2
epsilon = .1
psi_now = psi_hat(n = 500, pw, pz, full_model, reduced_model)


EIF(index, w[index], z[index], y[index], pw, pw_index,
   pz, pz_index, pwz, full_model, reduced_model, psi_now)
third_term_wz(index, w[index], z[index], y[index], pw, pw_index,
              pz, pz_index, pwz, full_model, reduced_model)




psi_curr_1 <- hat_psi(w, z, y, pw, pw_index, 
                      pz, pz_index, pwz, full_model, reduced_model, psi_now)
psi_curr_2_1 <- hat_psi(w, z, y, pw_2, pw_index_2_1, 
                        pz_2, pz_index_2_1, pwz_2_1, full_model_2, reduced_model_2, psi_now)
phi_curr = c(psi_curr_1, psi_curr_2_1 )

### Now, we can start by considering the one-step correction and its 
### corresponding coverage etc.
pw_index_update = pw_index
pz_index_update = pz_index
pwz_update = pwz




psi_now = psi_hat(n = 300, pw, pz, full_model, reduced_model)
tilde_phi <- c()
eif_vec = c()
decorr_vec = c()
TT1 = Sys.time()
for (index in 1:length(y)){
  psi_now = psi_hat(n = 300, pw, pz, full_model, reduced_model)
  tilde_phi[index] =
    (EIF(index, w[index], z[index], y[index], pw, pw_index,
         pz, pz_index, pwz, full_model, reduced_model, psi_now) * phi_curr[1] +
       EIF(index, w[index], z[index], y[index], pw_2, pw_index_2_1,
           pz_2, pz_index_2_1, pwz_2, full_model_2, reduced_model_2, psi_now) * phi_curr[2])/
    norm(phi_curr, "2")
  
  eif_vec[index] = EIF(index, w[index], z[index], y[index], pw_index_update, pw_index_update, 
                       pz_index_update, pz_index_update, pwz_update, full_model, reduced_model, psi_now)
  decorr_vec[index] = DecorrLOCO(index,  w[index], z[index], y[index], pw, pw_index, 
                                 pz, pz_index, pwz, full_model, reduced_model)

}
Sys.time()-TT1
eif_estimate = 1/2*mean(decorr_vec)
var_estimate = mean(eif_vec^2)


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

gamma = optimise(likelihood_nor, c(-1, 1))
epsilon = optimise(likelihood_ori, c(-1, 1))



### Next, we'll try to update our estimator, and see how the coverage looks like.
### To do this, we need to update the densities wrt the estimated values.
### A bit of reminder on what we need to update on:
### pw, pz, pw_index, pz_index, pwz

tilde_phi_applied = exp(tilde_phi*gamma$objective)
pw_index_update = pw_index
pw_index_update$estimate = (pw_index$estimate*tilde_phi_applied)/sum(pw_index$estimate*tilde_phi_applied)
pz_index_update = pz_index
pz_index_update$estimate = (pz_index$estimate*tilde_phi_applied)/sum(pz_index$estimate*tilde_phi_applied)
pwz_update = pwz
pwz_update$estimate = (pwz$estimate * tilde_phi_applied) / sum(pwz$estimate * tilde_phi_applied)



eif_update_vec = c()
decorr_update_vec = c()
for (index in 1:length(y)){
  eif_update_vec[index] = EIF(index, w[index], z[index], y[index], pw_index_update, pw_index_update, 
                              pz_index_update, pz_index_update, pwz_update, full_model, reduced_model, psi_hat(n = 200, pw_index_update, pz_index_update, full_model, reduced_model))
  decorr_update_vec[index] = DecorrLOCO(index,  w[index], z[index], y[index], pw_index_update, pw_index_update, 
                                        pz_index_update, pz_index_update, pwz_update, full_model, reduced_model)
}
1/2*mean(decorr_update_vec)

hat_phi_applied = exp(eif_vec*epsilon$objective)
pw_index_update_ori = pw_index
pw_index_update_ori$estimate = (pw_index$estimate*hat_phi_applied)/sum(pw_index$estimate*hat_phi_applied)
pz_index_update_ori = pz_index
pz_index_update_ori$estimate = (pz_index$estimate*hat_phi_applied)/sum(pz_index$estimate*hat_phi_applied)
pwz_update_ori = pwz
pwz_update_ori$estimate = (pwz$estimate * hat_phi_applied) / sum(pwz$estimate * hat_phi_applied)

eif_update_vec_ori = c()
decorr_update_vec_ori = c()
for (index in 1:length(y)){
  eif_update_vec_ori[index] = EIF(index, w[index], z[index], y[index], pw_index_update_ori, pw_index_update_ori, 
                                  pz_index_update_ori, pz_index_update_ori, pwz_update_ori, full_model, reduced_model, psi_hat(n = 100, pw_index_update_ori, pz_index_update_ori, full_model, reduced_model))
  decorr_update_vec_ori[index] = DecorrLOCO(index,  w[index], z[index], y[index], pw_index_update_ori, pw_index_update_ori, 
                                            pz_index_update_ori, pz_index_update_ori, pwz_update_ori, full_model, reduced_model)
}
1/2*mean(decorr_update_vec_ori)








sample_w = rnorm(100, 3, 1)
sample_z = rnorm(100, .9, 1.09)

# sample_w = sample(pw$eval.points, size = n, prob = pw$estimate)
# sample_z = sample(pz$eval.points, size = n, prob = pz$estimate)

predict_data_full = data.frame(w = sample_w, z = sample_z)
predict_data_reduced = data.frame(z = sample_z)
mean((predict(full_model, predict_data_full)-
        predict(reduced_model, predict_data_reduced))^2)

library(doParallel)
aaa = foreach(i =1:1000, .combine = "c")%do%{
  sample_w = rnorm(1000, 3, 1)
  sample_z = rnorm(1000, .9, 1.09)

  predict_data_full = data.frame(w = sample_w, z = sample_z)
  predict_data_reduced = data.frame(z = sample_z)
  mean((predict(full_model, predict_data_full)-
          predict(reduced_model, predict_data_reduced))^2)
}
mean(aaa)
