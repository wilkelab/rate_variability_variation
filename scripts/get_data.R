library(plyr)
library(dplyr)
library(tidyr)

enzyme_data = read.table(file = "../data_summaries/enzyme_dataset/jec_ddg_209_monomers.csv", sep = ",", header = T)
enzyme_out_data = data.frame(enzyme_data$pdb, enzyme_data$site, enzyme_data$rsa.tien, enzyme_data$r4s.jtt.ej)

names(enzyme_out_data) <- c("pdb_id", "site", "rsa", "r4s_JTT")

mono_data <- read.table("../data_summaries/enzyme_dataset/profiles_213_monomers.csv", sep  = ",", header = T)
enzyme_contacts <- mono_data %>% select(pdb, site, WCN, CN.13A.)
names(enzyme_contacts) <- c("pdb_id", "site", "cn", "wcn")
final_enzyme_data = inner_join(enzyme_contacts, enzyme_out_data)

final_enzyme_data = mutate(final_enzyme_data, data_id = "enzyme")
#write.csv(final_enzyme_data, "../data_summaries/output_files/enzyme_data.csv")

entropy_data = read.table("../data_summaries/enzyme_dataset/pdb_entropies.csv", sep = ",", header = T)

new_entropy_data = entropy_data %>% select(pdb, site, entropy_from_alignments)
names(new_entropy_data) <- c("pdb_id", "site", "entropy")
reduced_entropy_data = new_entropy_data %>% filter(pdb_id != "1BBS" | pdb_id != "1BS0" | pdb_id != "1DIN" | pdb_id != "1HPL")
#write.csv(reduced_entropy_data, "../data_summaries/output_files/enzyme_reduced_entropy_data.csv")

combo_enzyme_data = left_join(final_enzyme_data, reduced_entropy_data)
write.csv(combo_enzyme_data, "../data_summaries/output_files/all_enzyme_data.csv")


#Import Virus Data
hema_pdb = read.table("../data_summaries/virus_dataset/data_1RD8_AB.csv", sep = ",", header = T) #1RD8
hema_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/1RD8_AB_entropies.csv" , sep = ",", header = T)
hema_pdb  = hema_pdb %>% select(-entropy)
h_entropies =  hema_entropies %>% select(entropy)
hema_rates = read.table("../data_summaries/virus_dataset/mapped_rates/1RD8_AB_rates.csv" , sep = ",", header = T)
h_rates = hema_rates %>% select(rate)
hema_pdb = cbind(hema_pdb, h_entropies, h_rates)

dengue_pdb = read.table("../data_summaries/virus_dataset/data_2JLY_A.csv", sep = ",", header = T) #2JLY
dengue_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/2JLY_A_entropies.csv", sep = ",", header = T)
dengue_pdb  = dengue_pdb %>% select(-entropy)
d_entropies =  dengue_entropies %>% select(entropy)
dengue_rates = read.table("../data_summaries/virus_dataset/mapped_rates/2JLY_A_rates.csv" , sep = ",", header = T)
d_rates = dengue_rates %>% select(rate)
dengue_pdb = cbind(dengue_pdb, d_entropies, d_rates)


west_nile_pdb = read.table("../data_summaries/virus_dataset/data_2FP7_B.csv", sep = ",",  header = T) #2FP7
west_nile_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/2FP7_B_entropies.csv", sep = ",", header = T)
west_nile_pdb  = west_nile_pdb %>% select(-entropy)
wn_entropies =  west_nile_entropies %>% select(entropy)
west_nile_rates = read.table("../data_summaries/virus_dataset/mapped_rates/2FP7_B_rates.csv" , sep = ",", header = T)
wn_rates = west_nile_rates %>% select(rate)
west_nile_pdb = cbind(west_nile_pdb, wn_entropies, wn_rates)


jap_ence_pdb = read.table("../data_summaries/virus_dataset/data_2Z83_A.csv", sep = ",", header = T) #2Z83
jap_ence_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/2Z83_A_entropies.csv", sep = ",", header = T)
jap_ence_pdb  = jap_ence_pdb %>% select(-entropy)
je_entropies =  jap_ence_entropies %>% select(entropy)
je_rates = read.table("../data_summaries/virus_dataset/mapped_rates/2Z83_A_rates.csv" , sep = ",", header = T)
j_rates = je_rates %>% select(rate)
jap_ence_pdb = cbind(jap_ence_pdb, je_entropies, j_rates)


hep_c_pdb = read.table("../data_summaries/virus_dataset/data_3GOL_A.csv", sep = ",", header = T) #3GOL
hep_c_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/3GOL_A_entropies.csv", sep = ",", header = T)
hep_c_pdb  = hep_c_pdb %>% select(-entropy)
hc_entropies =  hep_c_entropies %>% select(entropy)
hep_c_rates = read.table("../data_summaries/virus_dataset/mapped_rates/3GOL_A_rates.csv" , sep = ",", header = T)
hc_rates = hep_c_rates %>% select(rate)
hep_c_pdb = cbind(hep_c_pdb, hc_entropies, hc_rates)


rift_valley_pdb = read.table("../data_summaries/virus_dataset/data_3LYF_A.csv", sep = ",",  header = T) #3LYF
rift_valley_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/3LYF_A_entropies.csv", sep = ",", header = T)
rift_valley_pdb  = rift_valley_pdb %>% select(-entropy)
rv_entropies =  rift_valley_entropies %>% select(entropy)
rift_valley_rates = read.table("../data_summaries/virus_dataset/mapped_rates/3LYF_A_rates.csv" , sep = ",", header = T)
rv_rates = rift_valley_rates %>% select(rate)
rift_valley_pdb = cbind(rift_valley_pdb, rv_entropies, rv_rates)

crim_congo_pdb = read.table("../data_summaries/virus_dataset/data_4AQF_B.csv",  sep = ",", header = T) #4AQF
crim_congo_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/4AQF_B_entropies.csv", sep = ",", header = T)
crim_congo_pdb = crim_congo_pdb %>% select(-entropy)
cc_entropies =  crim_congo_entropies %>% select(entropy)
crim_congo_rates = read.table("../data_summaries/virus_dataset/mapped_rates/4AQF_B_rates.csv" , sep = ",", header = T)
cc_rates = crim_congo_rates %>% select(rate)
crim_congo_pdb = cbind(crim_congo_pdb, cc_entropies, cc_rates)

marg_bd_pdb = read.table("../data_summaries/virus_dataset/data_4GHA_A.csv", sep = ",",  header = T) #4GHA
marg_bd_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/4GHA_A_entropies.csv", sep = ",", header = T)
marg_bd_pdb  = marg_bd_pdb %>% select(-entropy)
m_entropies =  marg_bd_entropies %>% select(entropy)
marg_rates = read.table("../data_summaries/virus_dataset/mapped_rates/4GHA_A_rates.csv" , sep = ",", header = T)
m_rates = marg_rates %>% select(rate)
marg_bd_pdb = cbind(marg_bd_pdb, m_entropies, m_rates)


flu_np_pdb = read.table("../data_summaries/virus_dataset/data_4IRY_A.csv", sep = ",", header = T ) # 4IRY
flu_np_entropies = read.table("../data_summaries/virus_dataset/mapped_entropies/4IRY_A_entropies.csv", sep = ",", header = T)
flu_np_pdb  = flu_np_pdb %>% select(-entropy)
flu_entropies =  flu_np_entropies %>% select(entropy)
flu_rates = read.table("../data_summaries/virus_dataset/mapped_rates/4IRY_A_rates.csv" , sep = ",", header = T)
f_rates = flu_rates %>% select(rate)
flu_np_pdb = cbind(flu_np_pdb, flu_entropies, f_rates)


virus_list = list(hema_pdb, dengue_pdb, west_nile_pdb, jap_ence_pdb, hep_c_pdb, rift_valley_pdb, crim_congo_pdb, marg_bd_pdb, flu_np_pdb)
virus_data = rbind(hema_pdb, dengue_pdb, west_nile_pdb, jap_ence_pdb, hep_c_pdb, rift_valley_pdb, crim_congo_pdb, marg_bd_pdb, flu_np_pdb)

virus_out_data <- virus_data %>% select(protein, res_num, rsa_cr, cn13_cr, wcn_cr, entropy, rate)
names(virus_out_data) <-c("pdb_id", "site", "rsa", "cn", "wcn", "entropy", "r4s_JTT")


virus_out_data <- mutate(virus_out_data, data_id = "virus")
write.csv(virus_out_data, "../data_summaries/output_files/all_virus_data.csv")

all_data <- rbind(combo_enzyme_data, virus_out_data) 
write.csv(all_data, "../data_summaries/output_files/all_pdb_data.csv")
