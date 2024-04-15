library(mvtnorm)
library(foreach)
n = 2*1e6
sig = diag(rep(1,10))


# nonlinear = foreach(rho = seq(.1, .9, .1), .combine = "c")%do%{
#   sig[1,2] = rho
#   sig[2,1] = rho
# 
#   X = rmvnorm(n, mean = rep(0,10), sigma = sig)
#   mean((lm(sin(X[,1])~X[,2])$residuals)^2)*25
# }

linear = 25*(1-seq(.1,.9,.1)^2)


linear_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-12/full_truth_linear_1000.csv")[, -1]




estimates_linear = abs(matrix(rowMeans(linear_data), byrow = F, nrow = 6)
                       - matrix(rep(linear,6), byrow = T, nrow= 6))
colnames(estimates_linear) <- paste0(seq(.1, .9, .1))
rownames(estimates_linear) <- c("tMLE", "CV-tMLE" ,"Plug_in", "Scaled", "naive-CV", "SS-tMLE")



library(ggplot2)
library(tidyr)

# Convert the matrix to a data frame
model_bias_df <- as.data.frame(estimates_linear)

# Add a column for the model numbers
model_bias_df$Model <- rownames(model_bias_df)

# Pivot longer so each row is a single observation
model_bias_long <- pivot_longer(model_bias_df, cols = -Model, names_to = "Rho", values_to = "Bias")

# Plot
ggplot(model_bias_long, aes(x = Rho, y = Bias, color = Model, group = Model)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Bias of tMLE vs Plugin under linear true model",
       x = "Rho",
       y = "Bias") + ylim(0,2)  
ggsave("./truth_linear.png", width = 1200, height = 936, units = "px", bg ="white")

