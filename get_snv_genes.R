library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

snv_frequency <- read_csv("data files/all_frequency.csv")
snv_frequency <- snv_frequency[complete.cases(snv_frequency), ]
snv_frequency$GBH_change <- snv_frequency$GBH_28 - snv_frequency$GBH_0
snv_frequency$CTRL_change <- snv_frequency$CTRL_28 - snv_frequency$CTRL_0
snv_frequency$GBH_CTRL_change <- snv_frequency$CTRL_change - snv_frequency$GBH_change
snv_frequency$GBH_CTRL_change_abs <- abs(snv_frequency$GBH_CTRL_change)
sig_changes <- subset(snv_frequency, GBH_CTRL_change_abs >= 0.7)

good_scaffolds <- read_csv("data files/good_mag_scaffolds.csv")

eggnog_genes <- read_tsv("data files/eggnog_genes.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("gene", "COG_ID", "Description", "Preferred_name")]
eggnog_genes$scaffold <- eggnog_genes$gene %>% str_extract("[^_]*_[^_]*")
background_cog <- subset(eggnog_genes, scaffold %in% good_scaffolds$scaffold)
background_cog <- subset(background_cog, is.na(COG_ID) == F)
#write.csv(background_cog, "data files/cog_background_genes.csv", row.names = F)

sig_snvs_sum <- snv_frequency %>% group_by(mag, gene) %>% summarize(snvs_in_gene = sum(GBH_CTRL_change_abs >= 0.7))
sig_snvs_sum <- subset(sig_snvs_sum, snvs_in_gene > 0)
sig_snvs_sum <- subset(sig_snvs_sum, !str_detect(gene, ","))

significant_genes <- left_join(sig_snvs_sum, background_cog) 
#write.csv(significant_genes, "data files/allele_shifts_significant_genes.csv", row.names = F)
significant_genes$one <- 1
significant_gene_sum <- significant_genes %>% group_by(mag) %>% summarize(sig_genes = sum(one))

bin_305_freq <- subset(sig_changes, mag == "bin.305")
bin_609_freq <- subset(sig_changes, mag == "bin.609")
bin_676_freq <- subset(sig_changes, mag == "bin.676")
