library(tidyverse)
library(dplyr)

setwd("/Users/Emma/Library/Mobile Documents/com~apple~CloudDocs/OneDrive Documents/phd docs/Chapter 1/Aim 1A")

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Jan30.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))
all_genes <- read_csv("MAG_gene_info_subsamp.csv")
all_genes <- subset(all_genes, select = -c(scaffold, timepoint, new_time, gene_length))
all_genes <- subset(all_genes, str_detect(new_name, "T2"))

for(MAG in mag_list){
  MAG_snvs <- subset(all_MAG_snvs, mag == MAG)
  MAG_genes <- subset(all_genes, mag == MAG)
  MAG_snvs_for_gene <- MAG_snvs[, c('new_name', 'gene', 'groups', 'final_ref_freq')] %>% subset(is.na(gene) == F)
  MAG_snvs_for_gene <- subset(MAG_snvs_for_gene, final_ref_freq < 1 & final_ref_freq > 0)
  MAG_snvs_for_gene$snv <- 1
  MAG_snvs_gene_sum <- MAG_snvs_for_gene %>% group_by(new_name, gene) %>% summarize(SNV_total = sum(snv))
  MAG_snvs_gene_sum <- subset(MAG_snvs_gene_sum, !str_detect(gene, ","))
  MAG_snvs_gene_sum <- subset(MAG_snvs_gene_sum, SNV_total < 300)
  
  MAG_genes <- subset(MAG_genes, coverage > 3)
  MAG_genes$rel_cov <- MAG_genes$coverage / MAG_genes$mag_coverage
  MAG_genes <- subset(MAG_genes, rel_cov < 3)
  MAG_genes <- subset(MAG_genes, select = c("gene", "coverage", "new_name"))
  MAG_genes_wide <- pivot_wider(MAG_genes, names_from = "new_name", values_from = "coverage")
  MAG_genes_wide <- MAG_genes_wide %>% na.omit()
  MAG_genes_long <- pivot_longer(MAG_genes_wide, cols = -c("gene"), names_to = "new_name", values_to = "coverage")
  
  MAG_gene_snvs <- left_join(MAG_genes_long[, c("gene", "new_name")], MAG_snvs_gene_sum[, c("gene", "new_name", "SNV_total")],  by = c("gene", "new_name"))
  MAG_gene_snvs$SNV_total <- with(MAG_gene_snvs, ifelse(is.na(MAG_gene_snvs$SNV_total), 0, SNV_total))
  gene_pop_matrix <- pivot_wider(MAG_gene_snvs, names_from = 'gene', values_from = "SNV_total")
  gene_pop_matrix$new_name <- gene_pop_matrix$new_name %>% str_sub(end = -7)
  gene_pop_matrix <- gene_pop_matrix[order(gene_pop_matrix$new_name),]
  write.csv(gene_pop_matrix, paste(MAG, "_gene_matrix_biallelic.csv", sep = ""), row.names = F)
  
}

for(MAG in mag_list){
  MAG_snvs <- subset(all_MAG_snvs, mag == MAG)
  MAG_genes <- subset(all_genes, mag == MAG)
  MAG_snvs_non_syn <- MAG_snvs[, c('new_name', 'gene', 'groups', 'final_ref_freq', 'mutation_type')] %>% subset(is.na(gene) == F)
  MAG_snvs_non_syn <- subset(MAG_snvs_non_syn, mutation_type == "N")
  MAG_snvs_non_syn$snv <- 1
  MAG_snvs_non_syn_sum <- MAG_snvs_non_syn %>% group_by(new_name, gene) %>% summarize(SNV_total = sum(snv))
  MAG_snvs_non_syn_sum <- subset(MAG_snvs_non_syn_sum, !str_detect(gene, ","))
  MAG_snvs_non_syn_sum <- subset(MAG_snvs_non_syn_sum, SNV_total < 300)
  
  MAG_genes <- subset(MAG_genes, coverage > 3)
  MAG_genes$rel_cov <- MAG_genes$coverage / MAG_genes$mag_coverage
  MAG_genes <- subset(MAG_genes, rel_cov < 3)
  MAG_genes <- subset(MAG_genes, select = c("gene", "coverage", "new_name"))
  MAG_genes_wide <- pivot_wider(MAG_genes, names_from = "new_name", values_from = "coverage")
  MAG_genes_wide <- MAG_genes_wide %>% na.omit()
  MAG_genes_long <- pivot_longer(MAG_genes_wide, cols = -c("gene"), names_to = "new_name", values_to = "coverage")
  
  MAG_gene_snvs_non_syn <- left_join(MAG_genes_long[, c("gene", "new_name")], MAG_snvs_non_syn_sum[, c("gene", "new_name", "SNV_total")],  by = c("gene", "new_name"))
  MAG_gene_snvs_non_syn$SNV_total <- with(MAG_gene_snvs_non_syn, ifelse(is.na(MAG_gene_snvs_non_syn$SNV_total), 0, SNV_total))
  gene_pop_matrix_non_syn <- pivot_wider(MAG_gene_snvs_non_syn, names_from = 'gene', values_from = "SNV_total")
  gene_pop_matrix_non_syn$new_name <- gene_pop_matrix_non_syn$new_name %>% str_sub(end = -7)
  gene_pop_matrix_non_syn <- gene_pop_matrix_non_syn[order(gene_pop_matrix_non_syn$new_name),]
  write.csv(gene_pop_matrix_non_syn, paste(MAG, "_gene_matrix_non_syn.csv", sep = ""), row.names = F)
  
}
