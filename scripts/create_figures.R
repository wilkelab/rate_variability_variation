library(dplyr)
library(ggplot2)

correlation_data = read.csv("../data/all_protein_summary_stats.csv")

p1 = ggplot(correlation_data, aes(x = mean_entropy, y = rsa_entropy_rho, color = data_id)) +
	theme_classic() + 
  	xlab("Mean Entropy") + 
  	ylab("Correlation RSA-Entropy") + 
  	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) 
  	
  	ggsave("../figures/mean_entropy_rsa_cor.pdf", plot  = p1, height = 4, width = 6)

p2 = ggplot(correlation_data, aes(x = var_entropy, y = rsa_entropy_rho, color = data_id)) +
	theme_classic() + 
  	xlab("Variance Entropy") + 
  	ylab("Correlation RSA-Entropy") + 
  	 theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) 
  	
  	ggsave("../figures/var_entropy_rsa_cor.pdf", plot = p2, height = 4, width = 6)
  

p3 = ggplot(correlation_data, aes(x = mean_entropy, y = cn_entropy_rho, color = data_id)) +
	theme_classic() + 
  	xlab("Mean Entropy") + 
  	ylab("Correlation CN-Entropy") + 
  	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) 
  	 
  	ggsave("../figures/mean_entropy_cn_cor.pdf", plot  = p3, height = 4, width = 6)
  	

p4 = ggplot(correlation_data, aes(x = var_entropy, y = cn_entropy_rho, color = data_id)) +
	theme_classic() + 
  	xlab("Variance Entropy") + 
  	ylab("Correlation CN-Entropy") + 
  	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) 
  	 
  	ggsave("../figures/var_entropy_cn_cor.pdf", plot  = p4, height = 4, width = 6)

p5 = ggplot(correlation_data, aes(x = mean_entropy, y = wcn_entropy_rho, color = data_id)) +
	theme_classic() + 
	xlab("Mean Entropy") + 
  	ylab("Correlation WCN-Entropy") +  
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) 
  	
	ggsave("../figures/mean_entropy_wcn_cor.pdf", plot = p5, height = 4, width = 6)

p6 = ggplot(correlation_data, aes(x = var_entropy, y = wcn_entropy_rho, color = data_id)) +
	theme_classic() +
	xlab("Variance Entropy") + 
  	ylab("Correlation WCN-Entropy") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.text.x = element_text(face="bold")) + 
  	theme(axis.line = element_line(size = 1.25)) + 
  	theme(axis.ticks = element_line(size = 1.25)) + 
  	geom_point(size = 2) 
  	
  	ggsave("../figures/var_entropy_wcn_cor.pdf", plot  = p6, height = 4, width = 6) 


  	
