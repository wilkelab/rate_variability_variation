library(dplyr)
HP_rates = read.table("../evol_rates/rate4site/rates/unnormalized/1RD8_AB_evol_rates.txt", sep = ",",  header = T) #1RD8_AB
HP_map = read.table("../pdb_maps/1RD8_AB_HP_map.txt", sep = "\t", header = T)
HP_data = cbind(HP_rates, HP_map) 
mod_HP_data = HP_data %>% filter(pdb_pos != "NA")
out_HP_data = mod_HP_data %>% select(site, rate)
write.csv(out_HP_data, "../data_summaries/virus_dataset/mapped_rates/1RD8_AB_rates.csv")

WNPB_rates = read.table("../evol_rates/rate4site/rates/unnormalized/2FP7_B_evol_rates.txt", sep = ",",  header = T) #2FP7_B
WNPB_map = read.table("../pdb_maps/2FP7_B_WNPB_map.txt", sep = "\t", header = T)
WNPB_data = cbind(WNPB_rates, WNPB_map) 
mod_WNPB_data = WNPB_data %>% filter(pdb_pos != "NA")
out_WNPB_data = mod_WNPB_data %>% select(site, rate)
write.csv(out_WNPB_data, "../data_summaries/virus_dataset/mapped_rates/2FP7_B_rates.csv")

CCHFN_rates = read.table("../evol_rates/rate4site/rates/unnormalized/4AQF_B_evol_rates.txt", sep = ",",  header = T) #4AQF_B
CCHFN_map = read.table("../pdb_maps/4AQF_B_CCHFN_map.txt", sep = "\t", header = T)
CCHFN_data = cbind(CCHFN_rates, CCHFN_map) 
mod_CCHFN_data = CCHFN_data %>% filter(pdb_pos != "NA")
out_CCHFN_data = mod_CCHFN_data %>% select(site, rate)
write.csv(out_CCHFN_data, "../data_summaries/virus_dataset/mapped_rates/4AQF_B_rates.csv")

DPH_rates = read.table("../evol_rates/rate4site/rates/unnormalized/2JLY_A_evol_rates.txt", sep = ",",  header = T) #2JLY_A
DPH_map = read.table("../pdb_maps/2JLY_A_DPH_map.txt", sep = "\t", header = T)
DPH_data = cbind(DPH_rates, DPH_map) 
mod_DPH_data = DPH_data %>% filter(pdb_pos != "NA")
out_DPH_data = mod_DPH_data %>% select(site, rate)
write.csv(out_DPH_data, "../data_summaries/virus_dataset/mapped_rates/2JLY_A_rates.csv")

HCP_rates = read.table("../evol_rates/rate4site/rates/unnormalized/3GOL_A_evol_rates.txt", sep = ",",  header = T) #3GOL_A
HCP_map = read.table("../pdb_maps/3GOL_A_HCP_map.txt", sep = "\t", header = T)
HCP_data = cbind(HCP_rates, HCP_map) 
mod_HCP_data = HCP_data %>% filter(pdb_pos != "NA")
out_HCP_data = mod_HCP_data %>% select(site, rate)
write.csv(out_HCP_data, "../data_summaries/virus_dataset/mapped_rates/3GOL_A_rates.csv")


INP_rates = read.table("../evol_rates/rate4site/rates/unnormalized/4IRY_A_evol_rates.txt", sep = ",",  header = T) #4IRY_A
INP_map = read.table("../pdb_maps/4IRY_A_INP_map.txt", sep = "\t", header = T)
INP_data = cbind(INP_rates, INP_map) 
mod_INP_data = INP_data %>% filter(pdb_pos != "NA")
out_INP_data = mod_INP_data %>% select(site, rate)
write.csv(out_INP_data, "../data_summaries/virus_dataset/mapped_rates/4IRY_A_rates.csv")


JEHN_rates = read.table("../evol_rates/rate4site/rates/unnormalized/2Z83_A_evol_rates.txt", sep = ",",  header = T) #2Z83_A
JEHN_map = read.table("../pdb_maps/2Z83_A_JEHN_map.txt", sep = "\t", header = T)
JEHN_data = cbind(JEHN_rates, JEHN_map) 
mod_JEHN_data = JEHN_data %>% filter(pdb_pos != "NA")
out_JEHN_data = mod_JEHN_data %>% select(site, rate)
write.csv(out_JEHN_data, "../data_summaries/virus_dataset/mapped_rates/2Z83_A_rates.csv")

MRNABD_rates = read.table("../evol_rates/rate4site/rates/unnormalized/4GHA_A_evol_rates.txt", sep = ",",  header = T) #4GHA_A
MRNABD_map = read.table("../pdb_maps/4GHA_A_MRNABD_map.txt", sep = "\t", header = T)
MRNABD_data = cbind(MRNABD_rates, MRNABD_map) 
mod_MRNABD_data = MRNABD_data %>% filter(pdb_pos != "NA")
out_MRNABD_data = mod_MRNABD_data %>% select(site, rate)
write.csv(out_MRNABD_data, "../data_summaries/virus_dataset/mapped_rates/4GHA_A_rates.csv")

RVFVNP_rates = read.table("../evol_rates/rate4site/rates/unnormalized/3LYF_A_evol_rates.txt", sep = ",",  header = T) #3LYF_A
RVFVNP_map = read.table("../pdb_maps/3LYF_A_RVFVNP_map.txt", sep = "\t", header = T)
RVFVNP_data = cbind(RVFVNP_rates, RVFVNP_map) 
mod_RVFVNP_data = RVFVNP_data %>% filter(pdb_pos != "NA")
out_RVFVNP_data = mod_RVFVNP_data %>% select(site, rate)
write.csv(out_RVFVNP_data, "../data_summaries/virus_dataset/mapped_rates/3LYF_A_rates.csv")

