library(ks)
library(np)
library(MASS)


### Set the random seed
set.seed(2010)


### We can start by including some helper functions
source("~/Working/Ning/Giles_Project_1/Code/Project Huitoukan/02-30/decorrEIF_setup_updated.R")


### Data generation: In principle, we should be using the Friedman equation 
### to estimate the stuff. 
rho = 0.3

w = rnorm(1000, 3, 1)
y = 2 *w + rnorm(length(w), 0, 1)
z = rho * w + rnorm(length(w), 0, 1)
x = cbind(w,z)
dist = (max(y) - min(y))
min = min(y)

y_s = y
y = (y - min(y))/(max(y) - min(y))


full_model = lm(y~1+w+z)
reduced_model = lm(y~1+z)

trunc_index = full_model$fitted.values>.005&full_model$fitted.values<(1-.005)
w = w[trunc_index]
y = y[trunc_index]
z = z[trunc_index]
x = x[trunc_index]
full_model$fitted.values = full_model$fitted.values[trunc_index]


data_wz = data.frame(w = w, z = z)



### Density estimates
pw_index = kde(w, eval.points = w, density = T)
pz_index = kde(z, eval.points = z, density = T)
pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)

### Conditional Density w|z
pwcz_bw = npcdensbw(data = data_wz, ydat = w, xdat = z)
pwcz = npcdens(bws = pwcz_bw, data = data_wz)

### Now, we are ready to conduct the update, notice that by construction,
### we only need to update on the tangent space of y|w,z and w|z.
psi_now = psi_hat(n = 200, pw_index, pz_index, full_model, reduced_model)


### And we'll start with the update on y|w,z
eif_h = c()

for (index in 1:length(y)){
  eif_h[index] = third_h(index, w[index], z[index],
                                   pw_index, pwcz, full_model, reduced_model)
}

regress_frame = data.frame(y = y, curr_offset = logit(full_model$fitted.values), obs = eif_h)
epsilon_model = glm(y ~ 0 + offset(curr_offset) + obs, data = regress_frame,  family = "quasibinomial")

full_model$fitted.values =  epsilon_model$fitted.values
update_full_model = full_model
update_full_model$fitted.values = epsilon_model$fitted.values 

### We'll now try to consider the updated version of the EIF

eif = c()
for (index in 1:length(y)){
  third_term <- third_update(index, w[index], z[index], y[index],
                              pw_index, pwcz, full_model, reduced_model,
                              update_full_model)
  first_term <- first_term_z(z[index], n = 100, pw_index, full_model, reduced_model)
  second_term <- second_term_w(w[index], n = 100, pz_index, full_model, reduced_model)
  eif[index] = 2*third_term + first_term + second_term
  
}

esitmate_now = mean(eif)/2*(dist^2)
var_now = mean((eif -2*psi_now)^2)*(dist^4)











