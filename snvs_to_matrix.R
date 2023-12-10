library(tidyverse)
library(dplyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Dec7.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))
all_genes <- read_csv("MAG_gene_info_subsamp.csv")
all_genes <- subset(all_genes, select = -c(scaffold, timepoint, new_time, gene_length))

for(MAG in mag_list){
  MAG_snvs <- subset(all_MAG_snvs, mag == MAG)
  MAG_genes <- subset(all_genes, mag == MAG)
  MAG_snvs_for_gene <- MAG_snvs[, c('name', 'gene', 'number_divergent', 'class', 'mutation_type')] %>% subset(is.na(gene) == F)
  MAG_snvs_for_gene <- subset(MAG_snvs_for_gene, class == "SNV")
  MAG_snvs_gene_sum <- MAG_snvs_for_gene %>% group_by(name, gene) %>% count(class = "SNV", name = "total")
  MAG_snvs_gene_sum <- subset(MAG_snvs_gene_sum, !str_detect(gene, ","))
  MAG_snvs_gene_sum <- subset(MAG_snvs_gene_sum, total < 300)
  
  MAG_genes <- subset(MAG_genes, coverage > 3)
  MAG_genes$rel_cov <- MAG_genes$coverage / MAG_genes$mag_coverage
  MAG_genes <- subset(MAG_genes, rel_cov < 3)
  MAG_genes <- subset(MAG_genes, select = c("gene", "coverage", "name"))
  MAG_genes_wide <- pivot_wider(MAG_genes, names_from = "name", values_from = "coverage")
  MAG_genes_wide <- MAG_genes_wide %>% na.omit()
  MAG_genes_long <- pivot_longer(MAG_genes_wide, cols = -c("gene"), names_to = "name", values_to = "coverage")
  
  MAG_gene_snvs <- left_join(MAG_genes_long[, c("gene", "name")], MAG_snvs_gene_sum[, c("gene", "name", "total")],  by = c("gene", "name"))
  MAG_gene_snvs$total <- with(MAG_gene_snvs, ifelse(is.na(MAG_gene_snvs$total), 0, total))
  gene_pop_matrix <- pivot_wider(MAG_gene_snvs, names_from = 'gene', values_from = "total")
  gene_pop_matrix <- gene_pop_matrix[order(gene_pop_matrix$name),]
  write.csv(gene_pop_matrix, paste(MAG, "_gene_no_SNS_matrix_subsamp.csv", sep = ""), row.names = F)
  
}
