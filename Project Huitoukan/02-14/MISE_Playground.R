library(ks)
set.seed(2024)
# Load necessary libraries
w = rnorm(1000, 3, 1)
y = 2 *w + rnorm(length(w), 0, 1)
z = rho * w + rnorm(length(w), 0, 1)
x = cbind(w,z)

pw_index = kde(w, eval.points = w, density = T)
density_estimate = kde(z, eval.points = z, density = T)


# Step 3: Define the True Density Function
true_density <- function(x) {
  dnorm(x, mean = 0.9, sd = 1.09)
}

# Step 4: Calculate the Integrated Squared Error (ISE) Numerically
# Interpolation of the KDE estimate to use in integration
kde_interpolate <- approxfun(density_estimate$eval.points, density_estimate$estimate)

# Function to compute squared error at each point
squared_error <- function(x) {
  (kde_interpolate(x) - true_density(x))^2
}

# Numerical integration of the squared error
ise_result <- integrate(squared_error, min(density_estimate$eval.points), max(density_estimate$eval.points))

# Print the ISE result
print(paste("ISE Value:", ise_result$value))

# Note: For MISE, you would repeat the KDE and ISE calculation for multiple bootstrap samples
# and average the ISE values. This step is not shown here due to its computational nature.

