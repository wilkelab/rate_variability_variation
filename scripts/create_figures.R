library(dplyr)
library(ggplot2)
library(cowplot)

correlation_data = read.csv("../data_summaries/output_files/all_protein_summary_stats.csv")
correlation_data = rename(correlation_data, dataset = data_id)
enzyme_data = filter(correlation_data, dataset == "enzyme")
selected_data = all_data %>% select(pdb_id, dataset, entropy) %>% filter(pdb_id == "1G24" | pdb_id == "1DDJ" | pdb_id == "1IU4" | pdb_id == "1HQC")
selected_data = rename(selected_data, PDB = pdb_id)

selected_cor_data = correlation_data %>% filter(pdb_id == "1G24" | pdb_id == "1DDJ" | pdb_id == "1IU4" | pdb_id == "1HQC")


all_data = read.csv("../data_summaries/output_files/all_pdb_data.csv")
all_data= rename(all_data, dataset = data_id)

p1 = ggplot(correlation_data, aes(x = mean_entropy, y = rsa_entropy_rho, color = dataset)) +
	theme_classic() + 
  	xlab("Mean Entropy") + 
  	ylab("Correlation RSA-Entropy") + 
  	theme(text = element_text(size = 18)) +  #Get rid of bold
  	theme(axis.line = element_line(size = 1.05)) + #Make this thinner
  	theme(axis.ticks = element_line(size = 1.05)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red"))
  	  	
  	ggsave("../figures/mean_entropy_rsa_cor.pdf", plot = p1, height = 4.5, width = 6)

p2 = ggplot(correlation_data, aes(x = var_entropy, y = rsa_entropy_rho, color = dataset)) +
	theme_classic() + 
  	xlab("Variance Entropy") + 
  	ylab("Correlation RSA-Entropy") + 
  	 theme(text = element_text(size = 18)) + 
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
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
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red"))   	 
  	ggsave("../figures/mean_entropy_cn_cor.pdf", plot  = p3, height = 4.5, width = 6)
  	

p4 = ggplot(correlation_data, aes(x = var_entropy, y = cn_entropy_rho, color = dataset)) +
	theme_classic() + 
  	xlab("Variance Entropy") + 
  	ylab("Correlation CN-Entropy") + 
  	theme(text = element_text(size = 18)) + 
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
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
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
  	geom_point(size = 2) +
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red")) 
  	
	ggsave("../figures/mean_entropy_wcn_cor.pdf", plot = p5, height = 4.5, width = 6)

p6 = ggplot(correlation_data, aes(x = var_entropy, y = wcn_entropy_rho, color = dataset)) +
	theme_classic() +
	xlab("Variance Entropy") + 
  	ylab("Correlation WCN-Entropy") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
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
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
  	geom_point(size = 2) +  
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
   	scale_color_manual(values = c("black", "red")) 

  	
  	ggsave("../figures/var_rate_rsa_cor.pdf", plot = p7, height = 4.5, width = 6) 

p8 = ggplot(correlation_data, aes(x = var_rate, y = rate_cn_rho, color = dataset)) +
	theme_classic() +
	xlab("Variance Rate4Site") + 
  	ylab("Correlation CN-Rate4Site") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
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
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
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
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
  	geom_point(size = 2) + 
  	stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	scale_color_manual(values = c("black", "red")) 
 	
  	ggsave("../figures/var_entropy_rate_cor.pdf", plot = p10, height = 4.5, width = 6) 
  	

p11 = ggplot(correlation_data, aes(x = mean_entropy, y = var_entropy, color = dataset)) +
	theme_classic() +
	xlab("Mean Entropy") + 
  	ylab("Variance Entropy") +
	theme(text = element_text(size = 18)) + 
  	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) + 
  	geom_point(data = selected_cor_data, aes(x = mean_entropy, y = var_entropy), shape = 1, size = 10) +
  	geom_point(data = enzyme_data, aes(x = mean_entropy, y = var_entropy), color = "black") + 
  	geom_point(size = 2) + 
  	#stat_smooth(method = "lm", se = F, data = enzyme_data) + 
  	annotate("text", x = 2.2, y = 0.3, label = "1G24", size = 4) + 
  	annotate("text", x = 1.75, y = 0.75, label = "1DD4", size = 4) + 
  	annotate("text", x = 0.45, y = 0.15, label = "1IU4", size = 4) +
  	annotate("text", x = 0.75, y = 0.6, label = "1HQC", size = 4) +
  	scale_color_manual(values = c("black", "red")) 
 	
 	
  	 ggsave("../figures/variance_entropy_plot.pdf", plot = p11, height = 4.5, width = 6) 
  	
p12 = ggplot(selected_data, aes(x = entropy, fill = PDB)) +
	geom_density(alpha  = 0.4) + 
	xlab("Entropy") + 
	ylab("Density") + 
	theme(text = element_text(size = 18)) +
	theme(axis.line = element_line(size = 1.05)) + 
  	theme(axis.ticks = element_line(size = 1.05)) +  
	scale_fill_manual(values = c("cyan", "blueviolet", "lightgreen", "deeppink")) 
	
	ggsave("../figures/density_plot.pdf", plot = p11, height = 4.5, width = 6) 

protein_plot = plot_grid(p11, p12, rows = 1, labels = c("A", "B"), label_size = 20)
ggsave("../figures/protein_ex_plot.pdf", plot = protein_plot, height = 4.5, width = 12)

	
	
#1G24 - (2.223767, 0.2334338)
#1DDJ - (1.499909, 0.7570813)
#1IU4 - (0.2270161, 0.1173657)
#1HQC - (0.7517428, 0.5424466)
	
