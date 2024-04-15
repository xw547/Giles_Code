### In this document, we'll try to consider the coverage of rf model 
### We'll consider two cases, where one of them would be
### 1) original method
### 2) enlarged CI as proposed by Wasserman et al.


### The reason why this is of interest was that that the performance of the 
### rf models under copula generating process shows superior results -- they 
### are 
### Hence, it would be meaningful to look at the coverage.

### We can start with considering the CV-tMLE, Scaled, and the naive_CV.

#source("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-05/copula_sample_testing.R")
linear_truth = 25*(1-seq(.1,.9,.1)^2)
truth_vec = rep(linear_truth, each = 3)
truth_mat = matrix(rep(truth_vec, each =200), nrow =27, byrow = T)


linear_data = read.csv("~/Working/Ning/Giles_Project_1/Code/LOCO CPI tMLE/2024-04-12/full_svm_linear_var_1000.csv")[, -1]


estimate_now = linear_data[sort(c(1 + 6*(0:8), 3 + 6*(0:8), 5 + 6*(0:8))), ]
var_estimate = linear_data[sort(c(2 + 6*(0:8), 4 + 6*(0:8), 6 + 6*(0:8))), ]


cover = rbind( truth_mat - sqrt(var_estimate/(2000-1))*qt(.975, 2000-1), truth_mat  + sqrt(var_estimate/(2000-1))*qt(.975, 2000-1))
CI_length = 2*rowMeans(sqrt(var_estimate/(2000-1))*qt(.95, 2000-1))

coverage = foreach(i = 1:27, .combine = c)%do%{
  upper_cover = cover[i+27,]
  lower_cover = cover[i,]
  mean(estimate_now[i,] < upper_cover & estimate_now[i,] > lower_cover)
}

coverage = matrix(coverage, nrow = 3)

library(ggplot2)
library(reshape2)


rownames(coverage) <- c("cv-TMLE", "scaled", "naiveCV")
colnames(coverage) <- seq(0.1, 0.9, by = 0.1)

# Convert the matrix to a data frame in long format
long_data <- melt(coverage, varnames = c("Method", "Sequence"), value.name = "Coverage")
long_data$Sequence <- as.numeric(as.character(long_data$Sequence))

# Plotting with ggplot2
ggplot(long_data, aes(x = Sequence, y = Coverage, group = Method, color = Method)) +
  geom_line() +
  geom_point() +
  labs(title = "Coverage of svm with copula", x = "Sequence", y = "Coverage") +
  theme_minimal() + ylim(0,1) + geom_hline(yintercept =  .95, color = "red")

ggsave("./svm_linear_coverage.png", width = 1200, height = 936, units = "px", bg ="white")


