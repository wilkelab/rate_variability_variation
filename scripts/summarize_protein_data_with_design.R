library(plyr)
library(dplyr)
enzyme_data  = read.csv("../data_summaries/output_files/all_enzyme_data.csv")
virus_data = read.csv("../data_summaries/output_files/all_virus_data.csv")

pdb_ids = levels(virus_data$pdb_id)
virus_summary = data.frame(pdb_id = character(0), cn_entropy_rho = numeric(0),cn_entropy_p_value = numeric(0), wcn_entropy_rho =numeric(0), wcn_entropy_p_value = numeric(0), rsa_entropy_rho = numeric(0),rsa_entropy_p_value = numeric(0), rate_wcn_rho = numeric(0), rate_wcn_p_value = numeric(0), rate_cn_rho = numeric(0), rate_cn_p_value = numeric(0), rate_rsa_rho = numeric(0), rate_rsa_p_value = numeric(0), mean_entropy = numeric(0), var_entropy = numeric(0), sd_entropy = numeric(0), var_rate = numeric(0), sd_rate = numeric(0), mean_rsa = numeric(0), mean_cn = numeric(0), mean_wcn = numeric(0), mean_designed_entropy = numeric(0), var_designed_entropy = numeric(0), sd_designed_entropy = numeric(0), design_entropy_rho = numeric(0), design_entropy_p_value = numeric(0), design_rate_rho = numeric(0), design_rate_p_value = numeric(0))

for(i in pdb_ids){

	data = filter(virus_data, pdb_id == i)
	r_cn = as.numeric(cor.test(data$entropy, data$cn, method = "spearman")$estimate)
	p_cn = cor.test(data$entropy, data$cn, method = "spearman")$p.value
	
	r_wcn = as.numeric(cor.test(data$entropy, data$wcn, method = "spearman")$estimate)
	p_wcn = cor.test(data$entropy, data$wcn, method = "spearman")$p.value
	
	r_rsa = as.numeric(cor.test(data$entropy, data$rsa, method = "spearman")$estimate)
	p_rsa = cor.test(data$entropy, data$rsa, method = "spearman")$p.value

	r2_cn = as.numeric(cor.test(data$r4s_JTT, data$cn, method = "spearman")$estimate)
	p2_cn = cor.test(data$r4s_JTT, data$cn, method = "spearman")$p.value
	
	r2_wcn = as.numeric(cor.test(data$r4s_JTT, data$wcn, method = "spearman")$estimate)
	p2_wcn = cor.test(data$r4s_JTT, data$wcn, method = "spearman")$p.value
	
	r2_rsa = as.numeric(cor.test(data$r4s_JTT, data$rsa, method = "spearman")$estimate)
	p2_rsa = cor.test(data$r4s_JTT, data$rsa, method = "spearman")$p.value
	
	r_design = as.numeric(cor.test(data$entropy, data$designed_entropy, method = "spearman")$estimate)
	p_design = cor.test(data$entropy, data$designed_entropy, method = "spearman")$p.value	

	r2_design = as.numeric(cor.test(data$r4s_JTT, data$cn, method = "spearman")$estimate)
	p2_design = cor.test(data$r4s_JTT, data$cn, method = "spearman")$p.value
	
	mean_ent = mean(data$entropy)
	var_ent = var(data$entropy)
	sd_ent = sd(data$entropy)
	
	var_evol_rate = var(data$r4s_JTT)
	sd_evol_rate = sd(data$r4s_JTT)
	
	mean_rsa = mean(data$rsa)
	mean_cn = mean(data$cn)
	mean_wcn = mean(data$wcn)
	
	mean_design_ent = mean(data$designed_entropy)
	var_design_ent = mean(data$designed_entropy)
	sd_design_ent = mean(data$designed_entropy)
	
	new_data = data.frame(pdb_id = i, cn_entropy_rho = r_cn, cn_entropy_p_value = p_cn, wcn_entropy_rho = r_wcn, wcn_entropy_p_value = p_wcn, rsa_entropy_rho = r_rsa, rsa_entropy_p_value = p_rsa, rate_wcn_rho = r2_wcn, rate_wcn_p_value = p2_wcn, rate_cn_rho = r2_cn, rate_cn_p_value = p2_cn, rate_rsa_rho = r2_rsa, rate_rsa_p_value = p2_rsa, mean_entropy = mean_ent, var_entropy = var_ent, sd_entropy = sd_ent, var_rate = var_evol_rate, sd_rate = sd_evol_rate, mean_rsa = mean_rsa, mean_cn = mean_cn, mean_wcn = mean_wcn, mean_designed_entropy = mean_design_ent, var_designed_entropy = var_design_ent, sd_designed_entropy = sd_design_ent, design_entropy_rho = r_design, design_entropy_p_value = p_design, design_rate_rho = r2_design, design_rate_p_value = p2_design)
	
	virus_summary = rbind(virus_summary, new_data)
}

reduced_enzyme_data = enzyme_data %>% filter(designed_entropy != "NA")
reduced_enzyme_data$pdb_id = factor(reduced_enzyme_data$pdb_id)
pdb_ids = levels(reduced_enzyme_data$pdb_id)
#pdb_ids = levels(enzyme_data$pdb_id)
enzyme_summary = data.frame(pdb_id = character(0), cn_entropy_rho = numeric(0),cn_entropy_p_value = numeric(0), wcn_entropy_rho = numeric(0), wcn_entropy_p_value = numeric(0), rsa_entropy_rho = numeric(0),rsa_entropy_p_value = numeric(0),
rate_wcn_rho = numeric(0), rate_wcn_p_value = numeric(0), rate_cn_rho = numeric(0), rate_cn_p_value = numeric(0), rate_rsa_rho = numeric(0), rate_rsa_p_value = numeric(0), mean_entropy = numeric(0), var_entropy = numeric(0), sd_entropy = numeric(0), var_rate = numeric(0), sd_rate = numeric(0), mean_rsa = numeric(0), mean_cn = numeric(0), mean_wcn = numeric(0),  mean_designed_entropy = numeric(0), var_designed_entropy = numeric(0), sd_designed_entropy = numeric(0), design_entropy_rho = numeric(0), design_entropy_p_value = numeric(0), design_rate_rho = numeric(0), design_rate_p_value = numeric(0))

for(i in pdb_ids){
	data = filter(reduced_enzyme_data, pdb_id == i)
	r_cn = as.numeric(cor.test(data$entropy, data$cn, method = "spearman")$estimate)
	p_cn = cor.test(data$entropy, data$cn, method = "spearman")$p.value
	
	r_wcn = as.numeric(cor.test(data$entropy, data$wcn, method = "spearman")$estimate)
	p_wcn = cor.test(data$entropy, data$wcn, method = "spearman")$p.value
	
	r_rsa = as.numeric(cor.test(data$entropy, data$rsa, method = "spearman")$estimate)
	p_rsa = cor.test(data$entropy, data$rsa, method = "spearman")$p.value

	r2_cn = as.numeric(cor.test(data$r4s_JTT, data$cn, method = "spearman")$estimate)
	p2_cn = cor.test(data$r4s_JTT, data$cn, method = "spearman")$p.value
	
	r2_wcn = as.numeric(cor.test(data$r4s_JTT, data$wcn, method = "spearman")$estimate)
	p2_wcn = cor.test(data$r4s_JTT, data$wcn, method = "spearman")$p.value
	
	r2_rsa = as.numeric(cor.test(data$r4s_JTT, data$rsa, method = "spearman")$estimate)
	p2_rsa = cor.test(data$r4s_JTT, data$rsa, method = "spearman")$p.value
	
	r_design = as.numeric(cor.test(data$entropy, data$designed_entropy, method = "spearman")$estimate)
	p_design = cor.test(data$entropy, data$designed_entropy, method = "spearman")$p.value	

	r2_design = as.numeric(cor.test(data$r4s_JTT, data$cn, method = "spearman")$estimate)
	p2_design = cor.test(data$r4s_JTT, data$cn, method = "spearman")$p.value
	
	mean_ent = mean(data$entropy)
	var_ent = var(data$entropy)
	sd_ent = sd(data$entropy)

	var_evol_rate = var(data$r4s_JTT)
	sd_evol_rate = sd(data$r4s_JTT)
	
	mean_rsa = mean(data$rsa)
	mean_cn = mean(data$cn)
	mean_wcn = mean(data$wcn)	

	mean_design_ent = mean(data$designed_entropy)
	var_design_ent = mean(data$designed_entropy)
	sd_design_ent = mean(data$designed_entropy)	
	
	new_data = data.frame(pdb_id = i, cn_entropy_rho = r_cn, cn_entropy_p_value = p_cn, wcn_entropy_rho = r_wcn, wcn_entropy_p_value = p_wcn, rsa_entropy_rho = r_rsa, rsa_entropy_p_value = p_rsa, rate_wcn_rho = r2_wcn, rate_wcn_p_value = p2_wcn, rate_cn_rho = r2_cn, rate_cn_p_value = p2_cn, rate_rsa_rho = r2_rsa, rate_rsa_p_value = p2_rsa, mean_entropy = mean_ent, var_entropy = var_ent, sd_entropy = sd_ent, var_rate = var_evol_rate, sd_rate = sd_evol_rate, mean_rsa = mean_rsa, mean_cn = mean_cn, mean_wcn = mean_wcn, mean_designed_entropy = mean_design_ent, var_designed_entropy = var_design_ent, sd_designed_entropy = sd_design_ent, design_entropy_rho = r_design, design_entropy_p_value = p_design, design_rate_rho = r2_design, design_rate_p_value = p2_design)
	
	enzyme_summary = rbind(enzyme_summary, new_data)
}

virus_summary = virus_summary %>% mutate(data_id = "virus")
enzyme_summary = enzyme_summary %>% mutate(data_id = "enzyme")
total_summary = rbind(enzyme_summary, virus_summary)

write.csv(total_summary, "../data_summaries/output_files/design_protein_summary_stats.csv")
