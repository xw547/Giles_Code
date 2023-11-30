library(ggplot2)
setwd("~/Working/Ning/Giles_Project_1/Code/Lu_Implementation/11-16")

rho_1 = read.table("./true_sim_160_withforest_1rho_200N.csv", sep=" ")
rho_1 = data.matrix(rho_1)

rho_7 = read.table("./true_sim_160_withforest_7rho_200N.csv", sep=" ")
rho_7 = data.matrix(rho_7)


coverage = rowMeans(rbind(rho_1[3, ], rho_7[3, ]))
mean_result = rowMeans(rbind(rho_1[1, ], rho_7[1, ]))
empirical_est = rowMeans(rbind(rho_1[2, ],rho_7[2, ]))
bias_result = mean_result - 2*(1-c(.1, .7)^2)
#plot(x = c(.1,.7), y = coverage)
#lines(x = c(.1,.7), y = coverage)
#plot(x = c(.1,.7), y = bias_result)
#lines(x = c(.1,.7), y = abs(bias_result))

plot_data = data.frame(rho = c(.1, .7), coverage = coverage, 
                       estimate = mean_result, not = 1-c(.1,.7)^2, 
                       twice = 2*(1-c(.1,.7)^2),
                       plus = (1-c(.1,.7)^2+.7))

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

