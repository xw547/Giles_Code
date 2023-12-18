### In this document, I'll present some results from the large sample illustration
### of exploration1, where we used lm to estimate \hat{y}

library("ggplot2")

setwd("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-15/")

explor2_rho1 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho1.csv"))[,-c(1,6,7)]
explor2_rho3 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho3.csv"))[,-c(1,6,7)]
explor2_rho5 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho5.csv"))[,-c(1,6,7)]
explor2_rho7 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho7.csv"))[,-c(1,6,7)]
explor2_rho9 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho9.csv"))[,-c(1,6,7)]

bias_estimation <- rbind(colMeans(explor2_rho1), colMeans(explor2_rho3),
                         colMeans(explor2_rho5), colMeans(explor2_rho7),
                         colMeans(explor2_rho9))
bias_estimation <- cbind(bias_estimation, rowSums(bias_estimation) - 2*bias_estimation[,3])
bias_estimation <- cbind(bias_estimation, 3-2*(seq(1,9,2)/10)^2, (seq(1,9,2)/10))
bias_estimation <- data.frame(bias_estimation)
colnames(bias_estimation) <- c("First_Term", "Second_Term", "Third_Term",
                               "Forth_Term", "Estimate", "True_Value", "rho")

# Define a custom color palette
my_colors <- c("#D16103", "#52854C", "#4E84C4", "#999999")
variable_names <- c("Second_Term", "Third_Term", "Forth_Term", "True_Value")
names(my_colors) <- variable_names

# Reshape the data into longer format
bias_estimation_long <- pivot_longer(bias_estimation, cols = c(Second_Term, Third_Term, Forth_Term, True_Value))

# Create the plot using the reshaped data and the custom color palette
estimation_plot <- ggplot(data = bias_estimation_long, aes(x = rho, y = value, color = name)) +
  geom_line() +
  scale_color_manual(values = my_colors) +
  ylab("Value") +
  labs(color = "Term")

estimation_plot = ggplot(data = bias_estimation, aes(x = rho)) + 
  geom_line(aes(x = rho, y = True_Value), color = "#999999") +
  geom_line(aes(x = rho, y = Estimate), color = "#CC79A7")


################################################################################
################################################################################
################################################################################
t_95 = qt(.975, 2000-1)
estimate_fun <- function(mat){
  return(rowSums(mat) - 2*mat[,3])
}


eif_mean = 3-2*seq(from =.1, to = .9, by = .2)^2

var_explor2_rho1 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho1.csv"))[, 6]
var_explor2_rho3 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho3.csv"))[, 6]
var_explor2_rho5 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho5.csv"))[, 6]
var_explor2_rho7 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho7.csv"))[, 6]
var_explor2_rho9 = data.matrix(read.csv("./crossfit_explore_2_2_both_rf_large_rho9.csv"))[, 6]

coverage_mat <- cbind(var_explor2_rho1, var_explor2_rho3, var_explor2_rho5, 
                      var_explor2_rho7, var_explor2_rho9)
estimate_mat <- cbind(estimate_fun(explor2_rho1), estimate_fun(explor2_rho3),
                      estimate_fun(explor2_rho5), estimate_fun(explor2_rho7),
                      estimate_fun(explor2_rho9))
SE_vec <- sqrt(colMeans(coverage_mat)/1999)

coverage <- c()

for(i in 1:5){
  eif_Lower = eif_mean[i] - t_95*SE_vec[i]/2
  eif_Upper = eif_mean[i] + t_95*SE_vec[i]/2
  eif_covered =  mean(estimate_mat[,i]<= eif_Upper&estimate_mat[,i]>=eif_Lower)
  coverage[i]  = eif_covered
}


coverage


bias_estimation <- cbind(bias_estimation, coverage)

coverage_plot = ggplot(data = bias_estimation, aes(x = rho)) + 
  geom_line(aes(x = rho, y = coverage), color = "#4E84C4") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") + ylim(0,1)

