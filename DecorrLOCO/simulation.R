library(foreach)
library(doParallel)


# Register parallel backend
# This example uses a local cluster, but you can adjust it for other backends like sockets or MPI
cl <- makeCluster(7)  # Here, we're using 2 cores
registerDoParallel(cl)

# Create a vector of numbers


# Use foreach to square each number in parallel
sim_result <- foreach(num = 1:196, .combine = "cbind") %dopar% {
  library(ks)
  #set.seed(2013)
  rho = 0.3
  
  w = rnorm(1000, 3, 1)
  y = 2 *w + rnorm(length(w), 0, 1)
  z = rho * w + rnorm(length(w), 0, 1)
  x = cbind(w,z)
  
  full_model = lm(y~1+w+z)
  reduced_model = lm(y~1+z)
  
  ### Density estimation:
  ### Notice that we are only interested in the density estimation of existing
  ### points.
  
  pw = density(w)
  pz = density(z)
  
  ### The following kde density estimates were used to estimate the density ratio 
  ### in the third term.
  
  pwz = kde(cbind(w,z), eval.points = cbind(w,z), density = T)
  pwz_grid = kde(cbind(w,z), density = T, gridsize = c(128L, 128L))
  argmin <- function(epsilon){
    return(likelihood(epsilon, w, z, y, pw, pw_index, 
                      pz, pz_index, pwz, pwz_grid, full_model, reduced_model, psi_now))
  }
  
  
  i = num
  epsilon_now = optimize(argmin, c(-0.5,0.5))$minimum
  c(DecorrLOCO(i, w[i], z[i], y[i], pw, pw_index, 
               pz, pz_index, pwz, full_model, reduced_model), 
    EIF(i, w[i], z[i], y[i], pw, a,  pz, b, pwz, full_model, reduced_model, 
        psi_hat(n = 100, pw, pz, full_model, reduced_model)),
    psi_hat(n = 100, pw, pz, full_model, reduced_model),epsilon_now )
  
  
}




 

