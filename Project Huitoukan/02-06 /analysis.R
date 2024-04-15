library(ggplot2)

data_PartOne <- read.csv("~/Working/Ning/Giles_Project_1/Code/Project Huitoukan/01-30/sim_rho3_updates.csv")
data_PartTwo <- read.csv("~/Working/Ning/Giles_Project_1/Code/Project Huitoukan/01-30/sim_rho3_updates_nor.csv")

data_Full <- na.omit(cbind(data_PartOne, data_PartTwo))[,-c(1,6)]
colnames(data_Full) <- c("1. One-Step", "2. Unweighted",
                         "3. One-Step Var", "4. Unweighted Var", "5. Weighted",
                         "6.Plug-in", "7. Weighted Var")

t_95 = qt(.975, 1000-1)

aaa = calculate_confidence_coverage(data_Full[,c(1,3)], 4.38)
bbb = calculate_confidence_coverage(data_Full[,c(2,4)], 4.38)
ccc = calculate_confidence_coverage(data_Full[,c(5,7)], 4.38)

calculate_confidence_coverage(data_Full[,c(1,3)], 4.38)
calculate_confidence_coverage(data_Full[,c(2,4)], 4.38)
calculate_confidence_coverage(data_Full[,c(5,7)], 4.38)

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
