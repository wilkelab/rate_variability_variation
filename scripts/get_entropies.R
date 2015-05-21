library(dplyr)
HP_entropies = read.table("../entropies/entropies/HP_entropies.txt", sep = "\t",  header = T) #1RD8_AB
HP_map = read.table("../entropies/pdb_maps/1RD8_AB_HP_map.txt", sep = "\t", header = T)
HP_data = cbind(HP_entropies, HP_map) 
mod_HP_data = HP_data %>% filter(pdb_pos != "NA")
out_HP_data = mod_HP_data %>% select(res_num, entropy)
write.csv(out_HP_data, "../data/mapped_entropies/1RD8_AB_entropies.csv")

WNPB_entropies = read.table("../entropies/entropies/WNPB_entropies.txt", sep = "\t",  header = T) #2FP7_B
WNPB_map = read.table("../entropies/pdb_maps/2FP7_B_WNPB_map.txt", sep = "\t", header = T)
WNPB_data = cbind(WNPB_entropies, WNPB_map) 
mod_WNPB_data = WNPB_data %>% filter(pdb_pos != "NA")
out_WNPB_data = mod_WNPB_data %>% select(res_num, entropy)
write.csv(out_WNPB_data, "../data/mapped_entropies/2FP7_B_entropies.csv")

CCHFN_entropies = read.table("../entropies/entropies/CCHFN_entropies.txt", sep = "\t",  header = T) #4AQF_B
CCHFN_map = read.table("../entropies/pdb_maps/4AQF_B_CCHFN_map.txt", sep = "\t", header = T)
CCHFN_data = cbind(CCHFN_entropies, CCHFN_map) 
mod_CCHFN_data = CCHFN_data %>% filter(pdb_pos != "NA")
out_CCHFN_data = mod_CCHFN_data %>% select(res_num, entropy)
write.csv(out_CCHFN_data, "../data/mapped_entropies/4AQF_B_entropies.csv")

DPH_entropies = read.table("../entropies/entropies/DPH_entropies.txt", sep = "\t",  header = T) #2JLY_A
DPH_map = read.table("../entropies/pdb_maps/2JLY_A_DPH_map.txt", sep = "\t", header = T)
DPH_data = cbind(DPH_entropies, DPH_map) 
mod_DPH_data = DPH_data %>% filter(pdb_pos != "NA")
out_DPH_data = mod_DPH_data %>% select(res_num, entropy)
write.csv(out_DPH_data, "../data/mapped_entropies/2JLY_A_entropies.csv")

HCP_entropies = read.table("../entropies/entropies/HCP_entropies.txt", sep = "\t",  header = T) #3GOL_A
HCP_map = read.table("../entropies/pdb_maps/3GOL_A_HCP_map.txt", sep = "\t", header = T)
HCP_data = cbind(HCP_entropies, HCP_map) 
mod_HCP_data = HCP_data %>% filter(pdb_pos != "NA")
out_HCP_data = mod_HCP_data %>% select(res_num, entropy)
write.csv(out_HCP_data, "../data/mapped_entropies/3GOL_A_entropies.csv")


INP_entropies = read.table("../entropies/entropies/INP_entropies.txt", sep = "\t",  header = T) #4IRY_A
INP_map = read.table("../entropies/pdb_maps/4IRY_A_INP_map.txt", sep = "\t", header = T)
INP_data = cbind(INP_entropies, INP_map) 
mod_INP_data = INP_data %>% filter(pdb_pos != "NA")
out_INP_data = mod_INP_data %>% select(res_num, entropy)
write.csv(out_INP_data, "../data/mapped_entropies/4IRY_A_entropies.csv")


JEHN_entropies = read.table("../entropies/entropies/JEHN_entropies.txt", sep = "\t",  header = T) #2Z83_A
JEHN_map = read.table("../entropies/pdb_maps/2Z83_A_JEHN_map.txt", sep = "\t", header = T)
JEHN_data = cbind(JEHN_entropies, JEHN_map) 
mod_JEHN_data = JEHN_data %>% filter(pdb_pos != "NA")
out_JEHN_data = mod_JEHN_data %>% select(res_num, entropy)
write.csv(out_JEHN_data, "../data/mapped_entropies/2Z83_A_entropies.csv")

MRNABD_entropies = read.table("../entropies/entropies/MRNABD_entropies.txt", sep = "\t",  header = T) #4GHA_A
MRNABD_map = read.table("../entropies/pdb_maps/4GHA_A_MRNABD_map.txt", sep = "\t", header = T)
MRNABD_data = cbind(MRNABD_entropies, MRNABD_map) 
mod_MRNABD_data = MRNABD_data %>% filter(pdb_pos != "NA")
out_MRNABD_data = mod_MRNABD_data %>% select(res_num, entropy)
write.csv(out_MRNABD_data, "../data/mapped_entropies/4GHA_A_entropies.csv")

RVFVNP_entropies = read.table("../entropies/entropies/RVFVNP_entropies.txt", sep = "\t",  header = T) #3LYF_A
RVFVNP_map = read.table("../entropies/pdb_maps/3LYF_A_RVFVNP_map.txt", sep = "\t", header = T)
RVFVNP_data = cbind(RVFVNP_entropies, RVFVNP_map) 
mod_RVFVNP_data = RVFVNP_data %>% filter(pdb_pos != "NA")
out_RVFVNP_data = mod_RVFVNP_data %>% select(res_num, entropy)
write.csv(out_RVFVNP_data, "../data/mapped_entropies/3LYF_A_entropies.csv")

