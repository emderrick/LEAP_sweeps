library(tidyverse)
library(dplyr)
library(cowplot)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs <- subset(all_MAG_snvs, new_time == 2)

gene_files <- list.files("95_profiles/",recursive = T, pattern=".*gene_info.tsv",full.names = T)

all_genes <- data.frame()
for(i in 1:length(gene_files)){
  pond_time_genes <- read.table(gene_files[i],sep="\t",header=T)
  timepoint <- gsub(".*profile_output/", "", gene_files[i]) %>% substr(1,9)
  pond_time_genes <- cbind(pond_time_genes,timepoint=rep(timepoint,nrow(pond_time_genes)))
  all_genes <- rbind(all_genes,pond_time_genes)
}

for(MAG in mag_list){
  MAG_snvs <- subset(all_MAG_snvs, mag == MAG)
  
  MAG_snvs_for_gene <- MAG_snvs[, c('name', 'gene', 'number_divergent', 'mutation_type')] %>% subset(is.na(gene) == F)
  MAG_snvs_for_gene$number_divergent <- with(MAG_snvs_for_gene, ifelse(is.na(MAG_snvs_for_gene$number_divergent), 0, number_divergent))
  MAG_snvs_gene_sum <- MAG_snvs_for_gene %>% group_by(name, gene) %>% summarize(total = sum(number_divergent))
  gene_pop_matrix <- pivot_wider(MAG_snvs_gene_sum, names_from = 'gene', values_from = "total")
  gene_pop_matrix[is.na(gene_pop_matrix)] = 0
  write.csv(gene_pop_matrix, paste(MAG, "_gene_matrix.csv", sep = ""), row.names = F)
  
  MAG_snvs_gene_nosyn <- subset(MAG_snvs_for_gene, mutation_type == "N")
  MAG_snvs_gene_nosyn_sum <- MAG_snvs_gene_nosyn %>% group_by(name, gene) %>% summarize(total = sum(number_divergent))
  gene_pop_matrix <- pivot_wider(MAG_snvs_gene_nosyn_sum, names_from = 'gene', values_from = "total")
  gene_pop_matrix[is.na(gene_pop_matrix)] = 0
  write.csv(gene_pop_matrix, paste(MAG, "_gene_nonsyn_matrix.csv", sep=""), row.names = F)
  
}
