library(plyr)
library(dplyr)
enzyme_data  = read.csv("../data/all_huang_et_al_data.csv")
virus_data = read.csv("../data/sharamoradi_et_al_data.csv")

pdb_ids = levels(virus_data$pdb_id)
virus_summary = data.frame(pdb_id = character(0), cn_entropy_rho =numeric(0),cn_entropy_p_value = numeric(0), wcn_entropy_rho =numeric(0), wcn_entropy_p_value = numeric(0), rsa_entropy_rho = numeric(0),rsa_entropy_p_value = numeric(0), mean_entropy = numeric(0), var_entropy = numeric(0), sd_entropy = numeric(0))

for(i in pdb_ids){

	data = filter(virus_data, pdb_id == i)
	r_cn = as.numeric(cor.test(data$entropy, data$cn, method = "spearman")$estimate)
	p_cn = cor.test(data$entropy, data$cn, method = "spearman")$p.value
	
	r_wcn = as.numeric(cor.test(data$entropy, data$wcn, method = "spearman")$estimate)
	p_wcn = cor.test(data$entropy, data$wcn, method = "spearman")$p.value
	
	r_rsa = as.numeric(cor.test(data$entropy, data$rsa, method = "spearman")$estimate)
	p_rsa = cor.test(data$entropy, data$rsa, method = "spearman")$p.value
	
	mean_ent = mean(data$entropy)
	var_ent = var(data$entropy)
	sd_ent = sd(data$entropy)
	
	new_data = data.frame(pdb_id = i, cn_entropy_rho = r_cn, cn_entropy_p_value = p_cn, wcn_entropy_rho = r_wcn, wcn_entropy_p_value = p_wcn, rsa_entropy_rho = r_rsa, rsa_entropy_p_value = p_rsa, mean_entropy = mean_ent, var_entropy = var_ent, sd_entropy = sd_ent)
	virus_summary = rbind(virus_summary, new_data)
}

pdb_ids = levels(enzyme_data$pdb_id)
enzyme_summary = data.frame(pdb_id = character(0), cn_entropy_rho =numeric(0),cn_entropy_p_value = numeric(0), wcn_entropy_rho = numeric(0), wcn_entropy_p_value = numeric(0), rsa_entropy_rho = numeric(0),rsa_entropy_p_value = numeric(0), mean_entropy = numeric(0), var_entropy = numeric(0), sd_entropy = numeric(0))

for(i in pdb_ids){
	data = filter(enzyme_data, pdb_id == i)
	r_cn = as.numeric(cor.test(data$entropy, data$cn, method = "spearman")$estimate)
	p_cn = cor.test(data$entropy, data$cn, method = "spearman")$p.value
	
	r_wcn = as.numeric(cor.test(data$entropy, data$wcn, method = "spearman")$estimate)
	p_wcn = cor.test(data$entropy, data$wcn, method = "spearman")$p.value
	
	r_rsa = as.numeric(cor.test(data$entropy, data$rsa, method = "spearman")$estimate)
	p_rsa = cor.test(data$entropy, data$rsa, method = "spearman")$p.value

	mean_ent = mean(data$entropy)
	var_ent = var(data$entropy)
	sd_ent = sd(data$entropy)
	
	new_data = data.frame(pdb_id = i, cn_entropy_rho = r_cn, cn_entropy_p_value = p_cn, wcn_entropy_rho = r_wcn, wcn_entropy_p_value = p_wcn, rsa_entropy_rho = r_rsa, rsa_entropy_p_value = p_rsa, mean_entropy = mean_ent, var_entropy = var_ent, sd_entropy = sd_ent)
	enzyme_summary = rbind(enzyme_summary, new_data)
}

virus_summary = virus_summary %>% mutate(data_id = "virus")
enzyme_summary = enzyme_summary %>% mutate(data_id = "enzyme")
total_summary = rbind(enzyme_summary, virus_summary)

write.csv(total_summary, "../data/all_protein_summary_stats.csv")


