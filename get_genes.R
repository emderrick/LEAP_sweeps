library(tidyverse)
library(dplyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

COG_gene_to_category <- read.table("cog-20.to_category.tsv", header = FALSE, sep = "\t")

eggnog_genes <- read_tsv("eggnog_output_oct/eggnog_genes_oct_nt.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$mag <- eggnog_genes$gene %>% substr(1,12)
eggnog_genes <- subset(eggnog_genes, mag %in% mag_list)
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("mag", "gene", "COG_ID", "Description", "Preferred_name")]
background_cog <- subset(eggnog_genes, is.na(COG_ID) == F)
write.csv(background_cog, "cog_background_genes.csv", row.names = F)

gene_changes_pass <- read_csv("gene_coverage_sig_genes_subsamp.csv")
gene_decrease <- subset(gene_changes_pass, cov_dif < -0.5)
gene_increase <- subset(gene_changes_pass, cov_dif > 0.5)
gene_cov_significant <- left_join(gene_changes_pass[, c("gene", "abs_val")], background_cog)
write.csv(gene_cov_significant, "gene_cov_significant_all_subsamp.csv", row.names = F)
gene_cov_sig_increase <- left_join(gene_increase[, c("gene", "abs_val")], background_cog)
write.csv(gene_cov_sig_increase, "gene_cov_sig_increase_all_subsamp.csv", row.names = F)
gene_cov_sig_decrease <- left_join(gene_decrease[, c("gene", "abs_val")], background_cog)
write.csv(gene_cov_sig_decrease, "gene_cov_sig_decrease_all_subsamp.csv", row.names = F)

not_strict_parevol_genes <- read_csv("loose_MAG_significant_genes_subsamp.csv")
not_strict_significant_genes <- left_join(not_strict_parevol_genes, background_cog)
write.csv(not_strict_significant_genes, "significant_genes_loose_all_subsamp.csv", row.names = F)

strict_parevol_genes <- read_csv("strict_MAG_significant_genes_subsamp.csv")
strict_significant_genes <- left_join(strict_parevol_genes, background_cog)
write.csv(strict_significant_genes, "significant_genes_strict_all_subsamp.csv", row.names = F)

threshold_snvs <- read_csv("threshold_snvs_subsamp.csv")
threshold_snvs <- subset(threshold_snvs, pass == "yes")
threshold_snvs_sum <- threshold_snvs %>% group_by(mag, scaffold, gene) %>% summarize(snvs_in_gene = sum(pass == "yes"))
threshold_significant_genes <- left_join(threshold_snvs_sum, background_cog) %>% subset(is.na(gene) == F)
write.csv(threshold_significant_genes, "threshold_significant_genes_all_subsamp.csv", row.names = F)

not_strict_parevol_sum <- not_strict_significant_genes %>% group_by(mag) %>% summarize(loose_pos_genes = sum(direction == "positive"), loose_neg_genes = sum(direction == "negative"))
strict_parevol_sum <- strict_significant_genes %>% group_by(mag) %>% summarize(strict_pos_genes = sum(direction == "positive"), strict_neg_genes = sum(direction == "negative"))
threshold_significant_genes$count <- 1
threshold_snvs_gene_sum  <- threshold_significant_genes  %>% group_by(mag) %>% summarize(allele_frequency_genes = sum(count))
sig_gene_summary <- right_join(strict_parevol_sum, not_strict_parevol_sum,  by = c("mag"))
sig_gene_summary <- right_join(sig_gene_summary, threshold_snvs_gene_sum, by = c("mag"))

MAG_overlap_loose <- data.frame()
for(MAG in mag_list){
  MAG_parevol_genes <- subset(not_strict_significant_genes, mag == MAG)
  MAG_threshold_genes <- subset(threshold_significant_genes, mag == MAG)
  MAG_parevol_gene_vector <- as.vector(MAG_parevol_genes$gene)
  MAG_threshold_gene_vector <- as.vector(MAG_threshold_genes$gene)
  overlaping_genes <- as.data.frame(intersect(MAG_parevol_genes$gene, MAG_threshold_genes$gene))
  overlaping_genes$mag <- MAG
  MAG_overlap_loose <- rbind(MAG_overlap_loose, overlaping_genes)
}

MAG_overlap_loose$count <- 1
overlap_summary_loose <- MAG_overlap_loose %>% group_by(mag) %>% summarise(loose_allele_frequency = sum(count))

MAG_overlap_strict <- data.frame()
for(MAG in mag_list){
  MAG_parevol_genes <- subset(strict_significant_genes, mag == MAG)
  MAG_threshold_genes <- subset(threshold_significant_genes, mag == MAG)
  MAG_parevol_gene_vector <- as.vector(MAG_parevol_genes$gene)
  MAG_threshold_gene_vector <- as.vector(MAG_threshold_genes$gene)
  overlaping_genes <- as.data.frame(intersect(MAG_parevol_genes$gene, MAG_threshold_genes$gene))
  MAG_overlap_strict <- rbind(MAG_overlap_strict, overlaping_genes)
}
MAG_overlap_strict$mag <- MAG_overlap_strict$`intersect(MAG_parevol_genes$gene, MAG_threshold_genes$gene)` %>% substr(1,12)
MAG_overlap_strict$count <- 1
overlap_summary_strict <- MAG_overlap_strict %>% group_by(mag) %>% summarise(strict_allele_frequency = sum(count))

sig_gene_summary <- left_join(sig_gene_summary, overlap_summary_loose, by = c("mag"))
sig_gene_summary <- left_join(sig_gene_summary, overlap_summary_strict, by = c("mag"))
sig_gene_summary[is.na(sig_gene_summary)] <- 0
write.csv(sig_gene_summary, "significant_genes_summary_subsamp.csv", row.names = F)

# strict_significant_genes_I4_MAG_00006 <- subset(strict_significant_genes, mag == "I4_MAG_00006")
# strict_significant_genes_I4_MAG_00006 <- left_join(strict_significant_genes_I4_MAG_00006, COG_gene_to_category, by = c("COG_ID" = "V1"))
# strict_significant_genes_I4_MAG_00006 <- rename(strict_significant_genes_I4_MAG_00006, COG_category = V2)
# write.csv(strict_significant_genes_I4_MAG_00006, "I4_MAG_00006_examples_strict.csv", row.names = F)
# 
# loose_significant_genes_I4_MAG_00006 <- subset(not_strict_significant_genes, mag == "I4_MAG_00006")
# loose_significant_genes_I4_MAG_00006 <- left_join(loose_significant_genes_I4_MAG_00006, COG_gene_to_category, by = c("COG_ID" = "V1"))
# loose_significant_genes_I4_MAG_00006 <- rename(loose_significant_genes_I4_MAG_00006, COG_category = V2)
# write.csv(loose_significant_genes_I4_MAG_00006, "I4_MAG_00006_examples_loose.csv", row.names = F)
# 
# threshold_significant_genes_I4_MAG_00006 <- subset(threshold_significant_genes, mag == "I4_MAG_00006")
# threshold_significant_genes_I4_MAG_00006 <- left_join(threshold_significant_genes_I4_MAG_00006, COG_gene_to_category, by = c("COG_ID" = "V1"))
# threshold_significant_genes_I4_MAG_00006 <- rename(threshold_significant_genes_I4_MAG_00006, COG_category = V2)
# write.csv(threshold_significant_genes_I4_MAG_00006, "I4_MAG_00006_examples_threshold.csv", row.names = F)
# 
# strict_significant_genes_L3_MAG_00058 <- subset(strict_significant_genes, mag == "L3_MAG_00058")
# strict_significant_genes_L3_MAG_00058 <- left_join(strict_significant_genes_L3_MAG_00058, COG_gene_to_category, by = c("COG_ID" = "V1"))
# strict_significant_genes_L3_MAG_00058 <- rename(strict_significant_genes_L3_MAG_00058, COG_category = V2)
# write.csv(strict_significant_genes_L3_MAG_00058, "L3_MAG_00058_examples_strict.csv", row.names = F)
# 
# loose_significant_genes_L3_MAG_00058 <- subset(not_strict_significant_genes, mag == "L3_MAG_00058")
# loose_significant_genes_L3_MAG_00058 <- left_join(loose_significant_genes_L3_MAG_00058, COG_gene_to_category, by = c("COG_ID" = "V1"))
# loose_significant_genes_L3_MAG_00058 <- rename(loose_significant_genes_L3_MAG_00058, COG_category = V2)
# write.csv(loose_significant_genes_L3_MAG_00058, "L3_MAG_00058_examples_loose.csv", row.names = F)
# 
# threshold_significant_genes_L3_MAG_00058 <- subset(threshold_significant_genes, mag == "L3_MAG_00058")
# threshold_significant_genes_L3_MAG_00058 <- left_join(threshold_significant_genes_L3_MAG_00058, COG_gene_to_category, by = c("COG_ID" = "V1"))
# threshold_significant_genes_L3_MAG_00058 <- rename(threshold_significant_genes_L3_MAG_00058, COG_category = V2)
# write.csv(threshold_significant_genes_L3_MAG_00058, "L3_MAG_00058_examples_threshold.csv", row.names = F)
# 
# loose_significant_genes_E <- left_join(not_strict_significant_genes, COG_gene_to_category, by = c("COG_ID" = "V1"))
# loose_significant_genes_E <- rename(loose_significant_genes_E, COG_category = V2)
# loose_significant_genes_E <- subset(loose_significant_genes_E, COG_category == "E")
# write.csv(loose_significant_genes_E, "loose_significant_genes_E.csv", row.names = F)
# 
# threshold_significant_genes_E <- left_join(threshold_significant_genes, COG_gene_to_category, by = c("COG_ID" = "V1"))
# threshold_significant_genes_E <- rename(threshold_significant_genes_E, COG_category = V2)
# threshold_significant_genes_E <- subset(threshold_significant_genes_E, COG_category == "E")
# write.csv(threshold_significant_genes_E, "threshold_significant_genes_E.csv", row.names = F)

# threshold_snvs_sum  <- na.omit(threshold_snvs_sum)
# threshold_genes <- as.vector(threshold_snvs_sum$gene)
# parevol_strict_genes <- as.vector(strict_parevol_genes$gene)
# parevol_not_strict_genes <- as.vector(not_strict_parevol_genes$gene)
# 
# all_sig_genes <- list(threshold_genes, parevol_strict_genes, parevol_not_strict_genes)
# 
# gene_VD <- ggVennDiagram(all_sig_genes, 
#                          category.names = c("Allele Freq", "Strict", "Loose"), 
#                          label = c("count"), edge_size = 0) +
#   scale_x_continuous(expand = expansion(mult = .2))+
#   scale_fill_gradient('Genes', high = "purple4", low = "thistle")
# 
# save_plot("gene_venndiagram.jpeg", gene_VD)