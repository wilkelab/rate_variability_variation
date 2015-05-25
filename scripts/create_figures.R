library(dplyr)
library(ggplot2)
library(cowplot)

correlation_data = read.csv("../data_summaries/output_files/all_protein_summary_stats.csv")
correlation_data = rename(correlation_data, dataset = data_id)
enzyme_data = filter(correlation_data, dataset == "enzyme")

p1 = ggplot(correlation_data, aes(x = mean_entropy, y = rsa_entropy_rho, color = dataset)) +
	theme_classic() + 
  	xlab("Mean Entropy") + 
  	ylab("Correlation RSA-Entropy") + 
  	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red"))
  	  	
  	ggsave("../figures/mean_entropy_rsa_cor.pdf", plot = p1, height = 4.5, width = 6)

p2 = ggplot(correlation_data, aes(x = var_entropy, y = rsa_entropy_rho, color = dataset)) +
	theme_classic() + 
  	xlab("Variance Entropy") + 
  	ylab("Correlation RSA-Entropy") + 
  	 theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) +
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red")) 
 
  	
	ggsave("../figures/var_entropy_rsa_cor.pdf", plot = p2, height = 4, width = 6)

  
mean_entropy_plot = plot_grid(p1, p2, rows = 1, labels = c("A", "B"),label_size = 20)
ggsave("../figures/entropy_rsa_cor.pdf", plot = mean_entropy_plot, height = 4.5, width = 12)

p3 = ggplot(correlation_data, aes(x = mean_entropy, y = cn_entropy_rho, color = dataset)) +
	theme_classic() + 
  	xlab("Mean Entropy") + 
  	ylab("Correlation CN-Entropy") + 
  	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red")) 

  	 
  	ggsave("../figures/mean_entropy_cn_cor.pdf", plot  = p3, height = 4.5, width = 6)
  	

p4 = ggplot(correlation_data, aes(x = var_entropy, y = cn_entropy_rho, color = dataset)) +
	theme_classic() + 
  	xlab("Variance Entropy") + 
  	ylab("Correlation CN-Entropy") + 
  	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) +
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
   	scale_color_manual(values = c("black", "red")) 
 	
  	
  	 
  	ggsave("../figures/var_entropy_cn_cor.pdf", plot  = p4, height = 4.5, width = 6)

mean_cn_plot = plot_grid(p3, p4, rows = 1, labels = c("A", "B"), label_size = 20)
ggsave("../figures/entropy_cn_cor.pdf", plot = mean_cn_plot, height = 4.5, width = 12)


p5 = ggplot(correlation_data, aes(x = mean_entropy, y = wcn_entropy_rho, color = dataset)) +
	theme_classic() + 
	xlab("Mean Entropy") + 
  	ylab("Correlation WCN-Entropy") +  
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) +
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red")) 

  	
  	
	ggsave("../figures/mean_entropy_wcn_cor.pdf", plot = p5, height = 4.5, width = 6)

p6 = ggplot(correlation_data, aes(x = var_entropy, y = wcn_entropy_rho, color = dataset)) +
	theme_classic() +
	xlab("Variance Entropy") + 
  	ylab("Correlation WCN-Entropy") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red")) 
  	
  	ggsave("../figures/var_entropy_wcn_cor.pdf", plot  = p6, height = 4.5, width = 6) 


mean_wcn_plot = plot_grid(p5, p6, rows = 1, labels = c("A", "B"), label_size = 20)
ggsave("../figures/entropy_wcn_cor.pdf", plot = mean_wcn_plot, height = 4.5, width = 12)


p7 = ggplot(correlation_data, aes(x = var_rate, y = rate_rsa_rho, color = dataset)) +
	theme_classic() +
	xlab("Variance Rate4Site") + 
  	ylab("Correlation RSA-Rate4Site") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) +  
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
   	scale_color_manual(values = c("black", "red")) 

  	
  	ggsave("../figures/var_rate_rsa_cor.pdf", plot = p7, height = 4.5, width = 6) 

p8 = ggplot(correlation_data, aes(x = var_rate, y = rate_cn_rho, color = dataset)) +
	theme_classic() +
	xlab("Variance Rate4Site") + 
  	ylab("Correlation CN-Rate4Site") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) +
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_y_continuous(limits = c(-0.8, 0.2)) + 
  	scale_color_manual(values = c("black", "red")) 

  	
  	ggsave("../figures/var_rate_cn_cor.pdf", plot = p8, height = 4.5, width = 6) 
	
p9 = ggplot(correlation_data, aes(x = var_rate, y = rate_wcn_rho, color = dataset)) +
	theme_classic() +
	xlab("Variance Rate4Site") + 
  	ylab("Correlation WCN-Rate4Site") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_y_continuous(limits = c(-0.8, 0.2)) + 
  	scale_color_manual(values = c("black", "red")) 
 	
  	ggsave("../figures/var_rate_wcn_cor.pdf", plot = p9, height = 4.5, width = 6) 
	
	
mean_rate_plot = plot_grid(p8, p9, rows = 1, labels = c("A", "B"), label_size = 20)
ggsave("../figures/rate_cor.pdf", plot = mean_rate_plot, height = 4.5, width = 12)

	
p10 = ggplot(correlation_data, aes(x = var_entropy, y = rate_wcn_rho, color = dataset)) +
	theme_classic() +
	xlab("Variance Entropy") + 
  	ylab("Correlation WCN-Rate4Site") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red")) 
 	
  	ggsave("../figures/var_entropy_rate_cor.pdf", plot = p10, height = 4.5, width = 6)  