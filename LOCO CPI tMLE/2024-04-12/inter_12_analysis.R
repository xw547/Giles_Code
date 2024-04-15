library(mvtnorm)
library(foreach)
n = 2*1e6
sig = diag(rep(1,10))




linear = 25*(1-seq(.1, .9, .1)^2)*2


svm_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-12/inter_12_svm_linear_1000.csv")[, -1]
lm_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-12/inter_12_lm_linear_1000.csv")[, -1]
rf_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-12/inter_12_rf_linear_1000.csv")[, -1]



estimates_svm = abs(matrix(rowMeans(svm_data), byrow = F, nrow = 6)
                    - matrix(rep(linear,6), byrow = T, nrow= 6))
colnames(estimates_svm) <- paste0(seq(.1, .9, .1))
rownames(estimates_svm) <- c("tMLE", "CV-tMLE" ,"Plug_in", "Scaled", "naive-CV", "SS-tMLE")



estimates_lm = abs(matrix(rowMeans(lm_data), byrow = F, nrow = 6)
                   - matrix(rep(linear,6), byrow = T, nrow= 6))
colnames(estimates_lm) <- paste0(seq(.1, .9, .1))
rownames(estimates_lm) <- c("tMLE", "CV-tMLE" ,"Plug_in", "Scaled", "naive-CV", "SS-tMLE")


estimates_rf = abs(matrix(rowMeans(rf_data), byrow = F, nrow = 6)
                   - matrix(rep(linear,6), byrow = T, nrow= 6))
colnames(estimates_rf) <- paste0(seq(.1, .9, .1))
rownames(estimates_rf) <- c("tMLE", "CV-tMLE" ,"Plug_in", "Scaled", "naive-CV", "SS-tMLE")




library(ggplot2)
library(tidyr)

# Convert the matrix to a data frame
model_bias_df <- as.data.frame(estimates_svm)

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
       y = "Bias")  
ggsave("./inter_12_svm_linear.png", width = 1200, height = 936, units = "px", bg ="white")






# Convert the matrix to a data frame
model_bias_df <- as.data.frame(estimates_lm)

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
       y = "Bias")  
ggsave("./inter_12_lm_linear.png", width = 1200, height = 936, units = "px", bg ="white")



# Convert the matrix to a data frame
model_bias_df <- as.data.frame(estimates_rf)

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
       y = "Bias")  
ggsave("./inter_12_rf_linear.png", width = 1200, height = 936, units = "px", bg ="white")
