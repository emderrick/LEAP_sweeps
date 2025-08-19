library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

snv_frequency <- read_csv("refined data files/refined_all_frequency.csv")
snv_frequency <- snv_frequency[complete.cases(snv_frequency), ]
snv_frequency$GBH_change <- snv_frequency$GBH_28 - snv_frequency$GBH_0
snv_frequency$CTRL_change <- snv_frequency$CTRL_28 - snv_frequency$CTRL_0
snv_frequency$GBH_CTRL_change <- snv_frequency$CTRL_change - snv_frequency$GBH_change
snv_frequency$GBH_CTRL_change_abs <- abs(snv_frequency$GBH_CTRL_change)
sig_changes <- subset(snv_frequency, GBH_CTRL_change_abs >= 0.7)
sig_nonsyn_changes <- subset(sig_changes, mutation_type == "N")
sig_syn_changes <- subset(sig_changes, mutation_type == "S")

good_scaffolds <- read_csv("refined data files/good_refined_mag_scaffolds.csv")

eggnog_genes <- read_tsv("refined data files/eggnog_genes.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("gene", "COG_ID", "Description", "Preferred_name")]
eggnog_genes$scaffold <- eggnog_genes$gene %>% str_extract("[^_]*_[^_]*")
background_cog <- subset(eggnog_genes, scaffold %in% good_scaffolds$scaffold)
background_cog <- subset(background_cog, is.na(COG_ID) == F)
write.csv(background_cog, "refined data files/refined_cog_background_genes.csv", row.names = F)

sig_snvs_sum <- sig_changes %>% group_by(mag, gene) %>% count()
sig_snvs_sum <- subset(sig_snvs_sum, !str_detect(gene, ","))

sig_nonsyn_sum <- sig_nonsyn_changes %>% group_by(mag, gene) %>% count()
sig_syn_sum <- sig_syn_changes %>% group_by(mag, gene) %>% count()

significant_genes <- left_join(sig_snvs_sum, background_cog) 
write.csv(significant_genes, "refined data files/refined_allele_shifts_significant_genes.csv", row.names = F)

significant_nonsyn_genes <- left_join(sig_nonsyn_sum, background_cog) 
write.csv(significant_nonsyn_genes, "refined data files/refined_nonsyn_allele_shifts_significant_genes.csv", row.names = F)

significant_syn_genes <- left_join(sig_syn_sum, background_cog) 
write.csv(significant_syn_genes, "refined data files/refined_syn_allele_shifts_significant_genes.csv", row.names = F)
