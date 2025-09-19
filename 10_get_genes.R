library(tidyverse)
library(patchwork)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

good_scaffolds <- read_csv("data files/good_refined_mag_scaffolds.csv")
eggnog_genes <- read_tsv("data files/eggnog_genes.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("gene", "COG_ID", "Description", "Preferred_name")]
eggnog_genes$scaffold <- eggnog_genes$gene %>% str_extract("[^_]*_[^_]*")
background_cog <- subset(eggnog_genes, scaffold %in% good_scaffolds$scaffold)
background_cog <- subset(background_cog, is.na(COG_ID) == F)
write.csv(background_cog, "data files/cog_background_genes.csv", row.names = F)

GBH_gene_loss <- read_csv("data files/GBH_gene_loss.csv")
GBH_gene_loss <- left_join(GBH_gene_loss, background_cog) 
write.csv(GBH_gene_loss, "data files/GBH_gene_loss_significant_genes.csv", row.names = F)

GBH_gene_gain <- read_csv("data files/GBH_gene_gain.csv")
GBH_gene_gain <- left_join(GBH_gene_gain, background_cog) 
write.csv(GBH_gene_gain, "data files/GBH_gene_gain_significant_genes.csv", row.names = F)

CTRL_gene_loss <- read_csv("data files/CTRL_gene_loss.csv")
CTRL_gene_loss <- left_join(CTRL_gene_loss, background_cog) 
write.csv(CTRL_gene_loss, "data files/CTRL_gene_loss_significant_genes.csv", row.names = F)

snv_frequency <- read_csv("data files/snv_frequency_changes.csv")

sig_b_changes <- subset(snv_frequency, case == "case_b")
sig_b_nonsyn_changes <- subset(sig_b_changes, mutation_type == "N")
sig_b_syn_changes <- subset(sig_b_changes, mutation_type == "S")

sig_b_snvs_sum <- sig_b_changes %>% group_by(mag, gene) %>% count()
sig_b_snvs_sum <- subset(sig_b_snvs_sum, !str_detect(gene, ","))

sig_b_nonsyn_sum <- sig_b_nonsyn_changes %>% group_by(mag, gene) %>% count()
sig_b_syn_sum <- sig_b_syn_changes %>% group_by(mag, gene) %>% count()

significant_b_genes <- left_join(sig_b_snvs_sum, background_cog) 
write.csv(significant_b_genes, "data files/allele_shifts_b_significant_genes.csv", row.names = F)

significant_b_nonsyn_genes <- left_join(sig_b_nonsyn_sum, background_cog) 
write.csv(significant_b_nonsyn_genes, "data files/nonsyn_allele_shifts_b_significant_genes.csv", row.names = F)

significant_b_syn_genes <- left_join(sig_b_syn_sum, background_cog) 
write.csv(significant_b_syn_genes, "data files/syn_allele_shifts_b_significant_genes.csv", row.names = F)

sig_a_changes <- subset(snv_frequency, case == "case_a")
sig_a_nonsyn_changes <- subset(sig_a_changes, mutation_type == "N")
sig_a_syn_changes <- subset(sig_a_changes, mutation_type == "S")

sig_a_snvs_sum <- sig_a_changes %>% group_by(mag, gene) %>% count()
sig_a_snvs_sum <- subset(sig_a_snvs_sum, !str_detect(gene, ","))

sig_a_nonsyn_sum <- sig_a_nonsyn_changes %>% group_by(mag, gene) %>% count()
sig_a_syn_sum <- sig_a_syn_changes %>% group_by(mag, gene) %>% count()

significant_a_genes <- left_join(sig_a_snvs_sum, background_cog) 
write.csv(significant_a_genes, "data files/allele_shifts_a_significant_genes.csv", row.names = F)

significant_a_nonsyn_genes <- left_join(sig_a_nonsyn_sum, background_cog) 
write.csv(significant_a_nonsyn_genes, "data files/nonsyn_allele_shifts_a_significant_genes.csv", row.names = F)

significant_a_syn_genes <- left_join(sig_a_syn_sum, background_cog) 
write.csv(significant_a_syn_genes, "data files/syn_allele_shifts_a_significant_genes.csv", row.names = F)



