library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(viridis)
library(cowplot)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
good_examples<- list("I4_MAG_00006", "L3_MAG_00058")

eggnog_genes <- read_tsv("eggnog_output_oct/eggnog_genes_oct_nt.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$mag <- eggnog_genes$gene %>% substr(1,12)
eggnog_genes <- subset(eggnog_genes, mag %in% mag_list)
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("mag", "gene", "COG_ID", "COG_category", "Description", "Preferred_name")]
background_cog <- subset(eggnog_genes, is.na(COG_ID) == F)
write.csv(background_cog, "cog_background_genes.csv", row.names = F)

not_strict_parevol_genes <- read_csv("loose_MAG_significant_genes.csv")
not_strict_significant_genes <- left_join(not_strict_parevol_genes, background_cog)
not_strict_significant_genes_good_ex <- subset(not_strict_significant_genes, mag %in% good_examples)
write.csv(not_strict_significant_genes, "significant_genes_loose.csv", row.names = F)

strict_parevol_genes <- read_csv("strict_MAG_significant_genes.csv")
strict_significant_genes <- left_join(strict_parevol_genes, background_cog)
strict_significant_genes_good_ex <- subset(strict_significant_genes, mag %in% good_examples)
write.csv(strict_significant_genes, "significant_genes_strict.csv", row.names = F)

threshold_snvs <- read_csv("threshold_snvs.csv")
threshold_snvs <- subset(threshold_snvs, pass=="yes")
threshold_snvs$snv_count <- 1
threshold_snvs_sum <- threshold_snvs %>% group_by(mag, scaffold, gene) %>% summarize(snvs_in_gene = sum(snv_count))
threshold_significant_genes <- left_join(threshold_snvs_sum, background_cog) %>% subset(is.na(gene) == F)
threshold_significant_genes_good_ex <- subset(threshold_significant_genes, mag %in% good_examples)
write.csv(threshold_significant_genes, "threshold_significant_genes.csv", row.names = F)

threshold_snvs_sum  <- na.omit(threshold_snvs_sum)
threshold_genes <- as.vector(threshold_snvs_sum$gene)
parevol_strict_genes <- as.vector(strict_parevol_genes$gene)
parevol_not_strict_genes <- as.vector(not_strict_parevol_genes$gene)

all_sig_genes <- list(threshold_genes, parevol_strict_genes, parevol_not_strict_genes)

gene_VD <- ggVennDiagram(all_sig_genes, 
              category.names = c("Allele Freq", "Strict", "Loose"), 
              label = c("count"), edge_size = 0) +
              scale_x_continuous(expand = expansion(mult = .2))+
              scale_fill_gradient('Genes', high = "purple4", low = "thistle")

save_plot("gene_venndiagram.jpeg", gene_VD)

not_strict_parevol_genes$count <- 1
not_strict_parevol_sum <- not_strict_parevol_genes %>% group_by(mag) %>% summarize(genes = sum(count))

strict_parevol_genes$count <- 1
strict_parevol_sum <- strict_parevol_genes %>% group_by(mag) %>% summarize(genes = sum(count))

threshold_snvs_sum$count <- 1
threshold_snvs_gene_sum  <- threshold_snvs_sum  %>% group_by(mag) %>% summarize(genes = sum(count))

mag_list2 <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                  "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00019")

MAG_overlap <- data.frame()
for(MAG in mag_list2){
  MAG_parevol_genes <- subset(not_strict_parevol_genes, mag == MAG)
  MAG_threshold_genes <- subset(threshold_snvs_sum, mag == MAG)
  MAG_parevol_gene_vector <- as.vector(MAG_parevol_genes$gene)
  MAG_threshold_gene_vector <- as.vector(MAG_threshold_genes$gene)
  overlaping_genes <- as.data.frame(intersect(MAG_parevol_genes$gene, MAG_threshold_genes$gene))
  overlaping_genes$mag <- MAG
  MAG_overlap <- rbind(MAG_overlap, overlaping_genes)
}

MAG_overlap$count <- 1
overlap_summary <- MAG_overlap %>% group_by(mag) %>% summarise(total_overlap = sum(count))

