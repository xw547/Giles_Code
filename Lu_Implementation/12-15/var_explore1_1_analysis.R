### In this document, I'll present some results from the large sample illustration
### of exploration1, where we used lm to estimate \hat{y}

library("ggplot2")

setwd("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/12-15/")

explor1_rho1 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho1.csv"))[,-c(1,6)]
explor1_rho3 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho3.csv"))[,-c(1,6)]
explor1_rho5 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho5.csv"))[,-c(1,6)]
explor1_rho7 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho7.csv"))[,-c(1,6)]
explor1_rho9 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho9.csv"))[,-c(1,6)]

bias_estimation <- rbind(colMeans(explor1_rho1), colMeans(explor1_rho3),
                         colMeans(explor1_rho5), colMeans(explor1_rho7),
                         colMeans(explor1_rho9))
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
t_95 = qt(.975, 2000-1)


eif_mean = 3-2*seq(from =.1, to = .9, by = .2)^2

var_explor1_rho1 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho1.csv"))[, 6]
var_explor1_rho3 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho3.csv"))[, 6]
var_explor1_rho5 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho5.csv"))[, 6]
var_explor1_rho7 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho7.csv"))[, 6]
var_explor1_rho9 = data.matrix(read.csv("./var_explore_1_1_lm_large_rho9.csv"))[, 6]

coverage_mat <- cbind(var_explor1_rho1, var_explor1_rho3, var_explor1_rho5, 
                      var_explor1_rho7, var_explor1_rho9)
estimate_mat <- cbind(estimate_fun(explor1_rho1), estimate_fun(explor1_rho3),
                      estimate_fun(explor1_rho5), estimate_fun(explor1_rho7),
                      estimate_fun(explor1_rho9))
SE_vec <- sqrt(colMeans(coverage_mat)/1999)

coverage <- c()

for(i in 1:5){
  eif_Lower = eif_mean[i] - t_95*SE_vec[i]
  eif_Upper = eif_mean[i] + t_95*SE_vec[i]
  eif_covered =  mean(estimate_mat[,i]<= eif_Upper&estimate_mat[,i]>=eif_Lower)
  coverage[i]  = eif_covered
}

coverage

estimate_fun <- function(mat){
  return(rowSums(mat) - 2*mat[,3])
}

bias_estimation <- cbind(bias_estimation, coverage)

coverage_plot = ggplot(data = bias_estimation, aes(x = rho)) + 
  geom_line(aes(x = rho, y = coverage), color = "#4E84C4") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") + ylim(0,1)



