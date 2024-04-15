library(ggplot2)

data_PartOne <- read.csv("~/Working/Ning/Giles_Project_1/Code/Project Huitoukan/02-20/sim_rho3_updates.csv")
data_PartTwo <- read.csv("~/Working/Ning/Giles_Project_1/Code/Project Huitoukan/02-20/sim_rho3_updates_nor.csv")

data_Full <- na.omit(cbind(data_PartOne, data_PartTwo))[,-c(1,6)]
colnames(data_Full) <- c("1. One-Step", "2. Unweighted",
                         "3. One-Step Var", "4. Unweighted Var", "5. Weighted",
                         "6.Plug-in", "7. Weighted Var")

t_95 = qt(.975, 1000-1)

aaa = calculate_confidence_coverage(data_Full[,c(1,3)], 4.36)
bbb = calculate_confidence_coverage(data_Full[,c(2,4)], 4.36)
ccc = calculate_confidence_coverage(data_Full[,c(5,7)], 4.36)

calculate_confidence_coverage(data_Full[,c(1,3)], 4.36)
calculate_confidence_coverage(data_Full[,c(2,4)], 4.36) 
calculate_confidence_coverage(data_Full[,c(5,7)], 4.36)


calculate_confidence_coverage <- function(data_matrix, true_value) {
  # Extract estimates and variances
  estimates <- data_matrix[, 1]
  variances <- data_matrix[, 2]
  
  # Calculate standard errors
  standard_errors <- sqrt(variances/(1000-1))
  
  # Z value for 95% confidence
  z_value <- qt(.975, 1000-1)
  
  # Calculate confidence intervals
  lower_bounds <- estimates - z_value * standard_errors
  upper_bounds <- estimates + z_value * standard_errors
  
  # Check if the true value is within the interval
  coverage_count <- sum(lower_bounds <= true_value & upper_bounds >= true_value)
  
  # Calculate coverage proportion
  coverage <- coverage_count / length(estimates)
  
  return(c(coverage, 2*mean(z_value * standard_errors)))
}



# library(doParallel)
# library(foreach)
# 
# coef = foreach(i =1:10000, .combine = "c")%do%{
#   set.seed(i + 320394)
#   w = rnorm(1000, 3, 1)
#   z =  rho * w + rnorm(length(w), 0, 1)
#   y = 2 *w + rnorm(length(w), 0, 1)
#   lm(y~0+z)$coefficients
# }
# mean(coef)



# Load necessary library
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
# Assuming this is your matrix (each row corresponds to one-step, unweighted, weighted)
# And each column corresponds to coverage and ci length respectively
matrix_data <- matrix(cbind(aaa, bbb, ccc), nrow = 3, byrow = TRUE,
                      dimnames = list(c("One-step", "Unweighted", "Weighted"),
                                      c("Coverage", "CI Length")))

df <- as.data.frame(matrix_data) %>%
  rownames_to_column("Method") %>%
  pivot_longer(-Method, names_to = "Metric", values_to = "Value")

# Convert 'Method' from a character to a factor to control order in the plot
df$Method <- factor(df$Method, levels = c("One-step", "Unweighted", "Weighted"))

# Plot for Coverage
ggplot(df %>% filter(Metric == "Coverage"), aes(x = Method, y = Value, group = 1)) +
  geom_line(color = "blue") +
  geom_point(size = 3, color = "blue") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "Coverage with true density implemented", x = "Method", y = "Coverage")

