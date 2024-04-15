library(MASS)
n = 2000000
sig = diag(rep(1,10))


# nonlinear = foreach(rho = seq(.1, .9, .1), .combine = "c")%do%{
#   sig[1,2] = rho
#   sig[2,1] = rho
# 
#   X = rmvnorm(n, mean = rep(0,10), sigma = sig)
#   mean((lm(sin(X[,1])~X[,2])$residuals)^2)*25
# }

copula_ff = c()
for (i in 1:9){
  rho =  0.1*i
  
  
  sig = diag(rep(1,10))
  sig[1,2] = rho
  sig[2,1] = rho
  n_1 = 1e6
  
  X = rmvcopula(n_1, 10, sig)
  X_reduced = X[,-1]
  beta = c(5, rep(0,9))
  Y = X%*%beta 
  
  
  full_data = data.frame(y = Y, X = X)
  colnames(full_data) <- c("output", paste0("X",1:10))
  reduced_data = data.frame(y = Y, X_reduced = X_reduced)
  colnames(reduced_data) <- c("output", paste0("X",2:10))
  
  
  X1.X.2_model = lm(X1~X2, data = full_data[1:n_1,])
  
  copula_ff[i] = 25*mean((full_data[1:n_1,2] -X1.X.2_model$fitted.values)^2)
  
}

linear = copula_ff  

linear_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-05/copula_lm_linear_1000.csv")[, -1]




estimates_linear = abs(rbind(matrix(rowMeans(linear_data), byrow = F, nrow = 4))
                       - matrix(rep(linear,4), byrow = T, nrow= 4))
colnames(estimates_linear) <- paste0(seq(.1, .9, .1))
rownames(estimates_linear) <- c("tMLE", "CV-tMLE" ,"Plug_in", "Scaled")


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
       y = "Bias")  
ggsave("./copula_lm_linear.png", width = 1200, height = 936, units = "px", bg ="white")


