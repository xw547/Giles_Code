library(ks)
library(np)

### Set the random seed
set.seed(2010)


### We can start by including some helper functions
source("~/Working/Ning/Giles_Project_1/Code/Project Huitoukan/02-14/decorrEIF_setup_updated.R")


### Data generation: In principle, we should be using the Friedman equation 
### to estimate the stuff. 
rho = 0.3

w = rnorm(1000, 3, 1)
y = 2 *w + rnorm(length(w), 0, 1)
z = rho * w + rnorm(length(w), 0, 1)
x = cbind(w,z)


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
data_ywz = data.frame(y = y, x = x)


### Density estimates
pw_index = kde(w, eval.points = w, density = T)
pz_index = kde(z, eval.points = z, density = T)
pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)

### Conditional Density w|z
pwcz_bw = npcdensbw(data = data_wz, ydat = w, xdat = z)
pwcz = npcdens(bws = pwcz_bw, data = data_wz)

### Conditional Density for y|w,z
pycwz_bw = npcdensbw(data = data_ywz, ydat = y, xdat = x)
pycwz = npcdens(bws = pycwz_bw, data = data_ywz)

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

update_full_model = full_model
update_full_model$fitted.values = epsilon_model$fitted.values


### Now we can proceed with the update on w|z
eif_second = c()
for (index in 1:length(y)){
  eif_second[index] = second_term_w(w[index], n = 100, pz_index, full_model, reduced_model)
}
epsilon_second = optimise(likelihood_first, c(-1/2, 1/2), y = y, pycwz = pwcz, eif_first = eif_second)

### z update: just a test
eif_third = c()
for (index in 1:length(y)){
  eif_third[index] = first_term_z(z[index], n = 100, pw_index, full_model, reduced_model)
}
epsilon_third = optimise(likelihood_ori, c(-1/2, 1/2), pwz = pz_index, eif_vec = eif_third, y = y)


### Updating the distributions, notice that in this case, we'll only need to update 
### p(w, z)
hat_phi_second = exp(eif_second * epsilon_second$objective)
pwcz$condens = (pwcz$condens * hat_phi_second)
pwz$estimate = pwcz$condens * pz_index$estimate
pz_index$estimate = pz_index$estimate *  exp(eif_third * epsilon_third$objective)














