library(tidyverse)
library(dplyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

eggnog_genes <- read_tsv("eggnog_output_oct/eggnog_genes_oct_nt.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$mag <- eggnog_genes$gene %>% substr(1,12)
eggnog_genes <- subset(eggnog_genes, mag %in% mag_list)
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("mag", "gene", "COG_ID", "Description", "Preferred_name")]
background_cog <- subset(eggnog_genes, is.na(COG_ID) == F)
write.csv(background_cog, "cog_background_genes.csv", row.names = F)

gene_changes_pass <- read_csv("gene_coverage_sig_genes_subsamp.csv")
gene_decrease <- subset(gene_changes_pass, T2_cov_dif < -0.5)
gene_increase <- subset(gene_changes_pass, T2_cov_dif > 0.5)
gene_cov_significant <- left_join(gene_changes_pass[, c("gene", "T2_cov_dif", "mag")], background_cog)
write.csv(gene_cov_significant, "gene_cov_significant_all_subsamp.csv", row.names = F)
gene_cov_sig_increase <- left_join(gene_increase[, c("gene", "T2_cov_dif", "mag")], background_cog)
write.csv(gene_cov_sig_increase, "gene_cov_sig_increase_all_subsamp.csv", row.names = F)
gene_cov_sig_decrease <- left_join(gene_decrease[, c("gene", "T2_cov_dif", "mag")], background_cog)
write.csv(gene_cov_sig_decrease, "gene_cov_sig_decrease_all_subsamp.csv", row.names = F)

not_strict_parevol_genes <- read_csv("loose_MAG_significant_genes_subsamp.csv")
not_strict_significant_genes <- left_join(not_strict_parevol_genes, background_cog)
write.csv(not_strict_significant_genes, "significant_genes_loose_all_subsamp.csv", row.names = F)

strict_parevol_genes <- read_csv("strict_MAG_significant_genes_subsamp.csv")
strict_significant_genes <- left_join(strict_parevol_genes, background_cog)
write.csv(strict_significant_genes, "significant_genes_strict_all_subsamp.csv", row.names = F)

threshold_snvs <- read_csv("threshold_snvs_subsamp.csv")
threshold_snvs <- subset(threshold_snvs, pass == "yes")
threshold_snvs_sum <- threshold_snvs %>% group_by(mag, scaffold, gene) %>% summarize(snvs_in_gene = sum(pass == "yes")) %>% subset(is.na(gene) == F)
threshold_snvs_sum <- subset(threshold_snvs_sum, !str_detect(gene, ","))
threshold_significant_genes <- left_join(threshold_snvs_sum, background_cog) 
write.csv(threshold_significant_genes, "threshold_significant_genes_all_subsamp.csv", row.names = F)

strict_significant_genes <- subset(strict_significant_genes, direction == "positive")
write.csv(strict_significant_genes, "parallel_decrease_genes.csv", row.names = F)

not_strict_significant_genes <- subset(not_strict_significant_genes, direction == "positive")
write.csv(not_strict_significant_genes, "SNV_decrease_genes.csv", row.names = F)

not_strict_parevol_sum <- not_strict_significant_genes %>% group_by(mag) %>% summarize(loose_pos_genes = sum(direction == "positive"))
strict_parevol_sum <- strict_significant_genes %>% group_by(mag) %>% summarize(strict_pos_genes = sum(direction == "positive"))

threshold_significant_genes$count <- 1
threshold_snvs_gene_sum  <- threshold_significant_genes  %>% group_by(mag) %>% summarize(allele_frequency_genes = sum(count))
coverage_sum <- read_csv("gene_cov_change_sum_subsamp.csv")
sig_gene_summary <- right_join(strict_parevol_sum, not_strict_parevol_sum,  by = c("mag"))
sig_gene_summary <- right_join(sig_gene_summary, threshold_snvs_gene_sum, by = c("mag"))

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
sig_gene_summary <- left_join(sig_gene_summary, overlap_summary_strict, by = c("mag"))
sig_gene_summary <- right_join(sig_gene_summary, coverage_sum)
sig_gene_summary[is.na(sig_gene_summary)] <- 0
write.csv(sig_gene_summary, "significant_genes_summary_subsamp.csv", row.names = F)


not_strict_parevol_NS_genes <- read_csv("loose_NS_sig_genes.csv")
not_strict_NS_genes <- left_join(not_strict_parevol_NS_genes, background_cog)
not_strict_NS_increase <- subset(not_strict_NS_genes, direction == "negative")
not_strict_NS_decrease <- subset(not_strict_NS_genes, direction == "positive")
write.csv(not_strict_NS_increase, "significant_genes_loose_inc_NS.csv", row.names = F)
write.csv(not_strict_NS_decrease, "significant_genes_loose_dec_NS.csv", row.names = F)

strict_parevol_NS_genes <- read_csv("strict_NS_sig_genes.csv")
strict_NS_genes <- left_join(strict_parevol_NS_genes, background_cog)
strict_NS_increase <- subset(strict_NS_genes, direction == "negative")
strict_NS_decrease <- subset(strict_NS_genes, direction == "positive")
write.csv(strict_NS_increase, "significant_genes_strict_inc_NS.csv", row.names = F)
write.csv(strict_NS_decrease, "significant_genes_strict_dec_NS.csv", row.names = F)


all_snv_info <- read_csv("all_MAG_SNVs_med_Jan30.csv")
all_snv_info <- subset(all_snv_info, select = c("groups", "new_name", "mutation_type"))
threshold_NS <- threshold_snvs[, c(4,3,7:19)]
threshold_NS <- pivot_longer(threshold_NS, cols = c(3:15), names_to = "new_name", values_to = "final_ref_freq", values_drop_na = T)
threshold_ns_type <- left_join(threshold_NS, all_snv_info)
threshold_ns_type <- threshold_ns_type[, c(1:3,5)]
threshold_NS_wide <- pivot_wider(threshold_ns_type, values_from = "mutation_type", names_from = "new_name")
