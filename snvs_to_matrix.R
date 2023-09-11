library(tidyverse)
library(dplyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs <- subset(all_MAG_snvs, new_time == 2)
MAG_snvs_wide <- read.csv("small_mag_hor_pass.csv")
snvs_pass <- subset(MAG_snvs_wide, pass == "yes")
snvs_pass$snv_GBH_count <- with(snvs_pass, ifelse(control == high, 1, 0))
snvs_pass$snv_CTRL_count <- with(snvs_pass, ifelse(control == low, 1, 0))
snvs_by_gene <- snvs_pass %>% group_by(name, gene) %>% summarize(total = sum(number_divergent))
  
for(MAG in mag_list){
  MAG_snvs <- subset(all_MAG_snvs, mag == MAG & number_divergent == 1 & is.na(gene) == F)
  MAG_snvs <- MAG_snvs[,c ('name', 'gene', 'number_divergent')]
  MAG_snvs_sum <- MAG_snvs %>% group_by(name, gene) %>% summarize(total = sum(number_divergent))
  gene_pop_matrix <- pivot_wider(MAG_snvs_sum, names_from = 'gene', values_from = "total")
  gene_pop_matrix[is.na(gene_pop_matrix)] = 0
  write.csv(gene_pop_matrix, paste(MAG, "_gene_pop_matrix.csv", sep=""), row.names = F)
}

for(MAG in mag_list){
  MAG_snvs <- subset(all_MAG_snvs, mag == MAG & number_divergent == 1 & is.na(gene) == F & mutation_type == "N")
  MAG_snvs <- MAG_snvs[,c ('name', 'gene', 'number_divergent')]
  MAG_snvs_sum <- MAG_snvs %>% group_by(name, gene) %>% summarize(total = sum(number_divergent))
  gene_pop_matrix <- pivot_wider(MAG_snvs_sum, names_from = 'gene', values_from = "total")
  gene_pop_matrix[is.na(gene_pop_matrix)] = 0
  write.csv(gene_pop_matrix, paste(MAG, "_gene_nonsyn_matrix.csv", sep=""), row.names = F)
}
