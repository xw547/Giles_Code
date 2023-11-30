library(ggplot2)
setwd("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/")

rho_1 = read.table("./sim_160_withforest_1_rho_100N.csv", sep=" ")
rho_1 = data.matrix(rho_1)

rho_3 <- read.table("./sim_160_withforest_3_rho_100N.csv", sep=" ")
rho_3 <- data.matrix(rho_3)

rho_5 <- read.table("./sim_160_withforest_5_rho_100N.csv", sep=" ")
rho_5 <- data.matrix(rho_5)

rho_7 <- read.table("./sim_160_withforest_7_rho_100N.csv", sep=" ")
rho_7 <- data.matrix(rho_7)

rho_9 <- read.table("./sim_160_withforest_9_rho_100N.csv", sep=" ")
rho_9 <- data.matrix(rho_9)


coverage = rowMeans(rbind(rho_1[3, ], rho_3[3, ], rho_5[3, ], rho_7[3, ], rho_9[3, ]))
mean_result = rowMeans(rbind(rho_1[1, ], rho_3[1, ], rho_5[1, ], rho_7[1, ], rho_9[1, ]))
empirical_est = rowMeans(rbind(rho_1[2, ], rho_3[2, ], rho_5[2, ], rho_7[2, ], rho_9[2, ]))
bias_result = mean_result - 2*(1-seq(.1, .9, by = .2)^2)
#plot(x = seq(.1, .9, by = .2), y = coverage)
#lines(x = seq(.1, .9, by = .2), y = coverage)
#plot(x = seq(.1, .9, by = .2), y = bias_result)
#lines(x = seq(.1, .9, by = .2), y = abs(bias_result))

plot_data = data.frame(rho = seq(.1, .9, by = .2), coverage = coverage, 
                       estimate = mean_result, not = 1-seq(.1, .9, by = .2)^2, 
                       twice = 2*(1-seq(.1, .9, by = .2)^2),
                       plus = (1-seq(.1, .9, by = .2)^2)+.75, 
                       ls = c(0.9791275, 0.8998923, 0.7413926, 0.5037403, 0.1872711),
                       forest = c(0.7559513, 0.7454708, 0.6614415, 0.4955970, 0.2403905)) 

coverage_plot = ggplot(data = plot_data, aes(x = rho) )+
  geom_line(aes(x = rho, y = coverage), color = "blue") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  xlab("rho") +
  ylab("coverage")+ 
  ggtitle("Coverage of 95% Confidence Interval") +
  theme(plot.title = element_text(hjust = 0.5))

estimate_plot = ggplot(data = plot_data, aes(x = rho))+  ylim(0,2)+
  geom_line(aes(x = rho, y = estimate, color = "estimate"), color = "blue") +
  geom_line(aes(x = rho, y = not, color = "1-rho^2"),  color = "red") +
  geom_line(aes(x = rho, y = twice, color = "2-2rho^2"),  color = "green") +
  geom_line(aes(x = rho, y = plus), color = "yellow")+
  xlab("rho") +
  ylab("Estimated CPI")+ 
  ggtitle("Average of EIF estimate with respect to rho") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color = "Legend Title")
   
mean_plot = ggplot(data = plot_data, aes(x = rho))+  ylim(0,2)+ 
  geom_line(aes(x = rho, y = not, color = "1-rho^2"),  color = "red")+
  geom_line(aes(x = rho, y = ls, color = "1-rho^2"),  color = "blue") +
  geom_line(aes(x = rho, y = forest, color = "1-rho^2"),  color = "green") +
  ggtitle("Plug in for E(f - f_-1) results") +
  ylab("Estimated Value")+
  theme(plot.title = element_text(hjust = 0.5))
  
