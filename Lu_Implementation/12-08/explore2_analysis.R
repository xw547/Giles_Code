### In this document, I'll present some results from the large sample illustration
### of exploration1, where we used lm to estimate \hat{y}

library("ggplot2")
library(matrixStats)
library(permimp)

setwd("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-08/")

explor2_rho1 = data.matrix(read.csv("./explore_1_rf_large_rho1.csv"))[,-1]
explor2_rho3 = data.matrix(read.csv("./explore_1_rf_large_rho3.csv"))[,-1]
explor2_rho5 = data.matrix(read.csv("./explore_1_rf_large_rho5.csv"))[,-1]
explor2_rho7 = data.matrix(read.csv("./explore_1_rf_large_rho7.csv"))[,-1]
explor2_rho9 = data.matrix(read.csv("./explore_1_rf_large_rho9.csv"))[,-1]

bias_estimation <- rbind(colMeans(explor2_rho1), colMeans(explor2_rho3),
                         colMeans(explor2_rho5), colMeans(explor2_rho7),
                         colMeans(explor2_rho9))
bias_estimation <- cbind(bias_estimation, rowSums(bias_estimation) - 2*bias_estimation[,3])
bias_estimation <- cbind(bias_estimation, 3-2*(seq(1,9,2)/10)^2, (seq(1,9,2)/10))
bias_estimation <- data.frame(bias_estimation)
colnames(bias_estimation) <- c("First_Term", "Second_Term", "Third_Term",
                               "Forth_Term", "Estimate", "True_Value", "rho")

estimation_plot = ggplot(data = bias_estimation, aes(x = rho))   +
  geom_line(aes(x = rho, y = Second_Term), color = "#D16103") +
  geom_line(aes(x = rho, y = Third_Term), color = "#52854C")  + 
  geom_line(aes(x = rho, y = Forth_Term), color = "#4E84C4")  +
  geom_line(aes(x = rho, y = True_Value), color = "#999999") 

estimation_plot = ggplot(data = bias_estimation, aes(x = rho)) + 
  geom_line(aes(x = rho, y = True_Value), color = "#999999") +
  geom_line(aes(x = rho, y = Estimate), color = "#CC79A7")



################################################################################
################################################################################
################################################################################
estimate_matrix <- t(rbind(estimate_fun(explor2_rho1), estimate_fun(explor2_rho3),
                         estimate_fun(explor2_rho5), estimate_fun(explor2_rho7),
                         estimate_fun(explor2_rho9)))

SE <- colMeans(estimate_matrix^2- colMeans(estimate_matrix)^2)

estimate_fun <- function(mat){
  return(rowSums(mat) - 2*mat[,3])
}

