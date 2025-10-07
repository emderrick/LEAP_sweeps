library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00110_1", "MAG_00179_1", "MAG_00194_1", "MAG_00197_1", "MAG_00201_1", "MAG_00674_1")

mag_genes <- read_csv("data files/MAG_gene_info.csv")

mag_genes$rel_cov <- mag_genes$gene_coverage / mag_genes$mag_coverage
mag_genes_wide <- pivot_wider(mag_genes[, c("mag", "gene", "Name_Time", "rel_cov")], names_from = "Name_Time", values_from = "rel_cov")
mag_genes_wide <- mag_genes_wide[-which(mag_genes_wide[3:18] > 3, arr.ind = TRUE)[, 1], ]

mag_genes_wide$CTRL_0 <-rowMeans(mag_genes_wide[, grep("Day 0 CTRL", colnames(mag_genes_wide))], na.rm = T)
mag_genes_wide$GBH_0 <- rowMeans(mag_genes_wide[, grep("Day 0 GBH", colnames(mag_genes_wide))], na.rm = T)
mag_genes_wide$CTRL_28 <- rowMeans(mag_genes_wide[, grep("Day 28 CTRL", colnames(mag_genes_wide))], na.rm = T)
mag_genes_wide$GBH_28 <- rowMeans(mag_genes_wide[, grep("Day 28 GBH", colnames(mag_genes_wide))], na.rm = T)
mag_gene_change <- mag_genes_wide[, c(1,2,19:22)]
mag_gene_change <- mag_gene_change %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))

mag_gene_change <- mag_gene_change[complete.cases(mag_gene_change), ]
mag_gene_change$GBH_change <- mag_gene_change$GBH_28 - mag_gene_change$GBH_0
mag_gene_change$CTRL_change <- mag_gene_change$CTRL_28 - mag_gene_change$CTRL_0
mag_gene_change$GBH_CTRL_change <- mag_gene_change$GBH_change - mag_gene_change$CTRL_change

mag_gene_change$GBH_change_abs <- abs(mag_gene_change$GBH_change)
mag_gene_change$CTRL_change_abs <- abs(mag_gene_change$CTRL_change)
mag_gene_change$GBH_CTRL_change_abs <- abs(mag_gene_change$GBH_CTRL_change)

gene_loss_GBH <- subset(mag_gene_change, CTRL_0 >= 0.6 & GBH_0 >= 0.6 & GBH_28 <= 0.1 & CTRL_28 >= 0.6)
write.csv(gene_loss_GBH, "data files/GBH_gene_loss.csv", row.names = F)

gene_gain_GBH <- subset(mag_gene_change, GBH_change >= 0.5 & CTRL_change < 0)
write.csv(gene_gain_GBH, "data files/GBH_gene_gain.csv", row.names = F)

gene_loss_CTRL <- subset(mag_gene_change, CTRL_0 >= 0.6 & GBH_0 >= 0.6 & GBH_28 >= 0.6 & CTRL_28 <= 0.1)
write.csv(gene_loss_CTRL, "data files/CTRL_gene_loss.csv", row.names = F)

eggnog_genes <- read_tsv("data files/eggnog_genes.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("gene", "COG_ID", "Description", "Preferred_name")]
eggnog_genes$scaffold <- eggnog_genes$gene %>% str_extract("[^_]*_[^_]*")
background_cog <- subset(eggnog_genes, gene %in% mag_genes_wide$gene)
background_cog <- subset(background_cog, is.na(COG_ID) == F)
write.csv(background_cog, "data files/cog_background_copy_number_genes.csv", row.names = F)

gene_loss_GBH <- left_join(gene_loss_GBH, background_cog) 
write.csv(gene_loss_GBH, "data files/GBH_gene_loss_significant_genes.csv", row.names = F)

gene_gain_GBH <- left_join(gene_gain_GBH, background_cog) 
write.csv(gene_gain_GBH, "data files/GBH_gene_gain_significant_genes.csv", row.names = F)

gene_loss_CTRL <- left_join(gene_loss_CTRL, background_cog) 
write.csv(gene_loss_CTRL, "data files/CTRL_gene_loss_significant_genes.csv", row.names = F)

sum_gain <- gene_gain_GBH %>% group_by(mag) %>% count()
write.csv(sum_gain, "data files/sum_GBH_gain.csv", row.names = F)

sum_loss <- gene_loss_GBH %>% group_by(mag) %>% count()
write.csv(sum_loss, "data files/sum_GBH_loss.csv", row.names = F)

