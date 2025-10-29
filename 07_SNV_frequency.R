library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

all_mag_SNVs <- read_csv("data files/MAG_SNV_depth_info.csv")
all_mag_SNVs <- subset(all_mag_SNVs, mag_coverage >= 5 & mag_breadth >= 0.7)

all_mag_SNVs$time <- with(all_mag_SNVs, ifelse(time == "1", "Day 0", "Day 28"))
all_mag_SNVs$pond_time <- paste(all_mag_SNVs$time, all_mag_SNVs$Name, sep = " ")
mag_snvs <- all_mag_SNVs[, c("gene", "pond_time", "mag", "group", "mutation_type", "final_ref_freq")]
postions_keep <- mag_snvs %>% group_by(group) %>% count(mutation_type)
postions_keep <- spread(postions_keep, mutation_type, n)
postions_keep$na_count <- rowSums(is.na(postions_keep[, c("I", "M", "N", "S")]))
postions_keep <- subset(postions_keep, na_count >= 3)
mag_snvs <- subset(mag_snvs, group %in% postions_keep$group)
mag_snvs <- mag_snvs %>% group_by(group) %>% fill(mutation_type, .direction = "updown") %>% ungroup()
mag_snvs_wide <- spread(mag_snvs, pond_time, final_ref_freq)

mag_snvs_wide$CTRL_0 <-rowMeans(mag_snvs_wide[, grep("Day 0 CTRL", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide$GBH_0 <- rowMeans(mag_snvs_wide[, grep("Day 0 GBH", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide$CTRL_28 <- rowMeans(mag_snvs_wide[, grep("Day 28 CTRL", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide$GBH_28 <- rowMeans(mag_snvs_wide[, grep("Day 28 GBH", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide <- mag_snvs_wide %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
mag_snvs_wide <- mag_snvs_wide[, c("mag", "gene", "group", "mutation_type", "CTRL_0", "GBH_0", "CTRL_28", "GBH_28")]
#write.csv(mag_snvs_wide, "data files/all_snv_frequency.csv", row.names = F)

snv_frequency <- mag_snvs_wide[complete.cases(mag_snvs_wide), ]
snv_frequency$GBH_change <- snv_frequency$GBH_28 - snv_frequency$GBH_0
snv_frequency$CTRL_change <- snv_frequency$CTRL_28 - snv_frequency$CTRL_0
snv_frequency$GBH_CTRL_change <- snv_frequency$GBH_change - snv_frequency$CTRL_change
snv_frequency$GBH_CTRL_T1 <- snv_frequency$GBH_0 - snv_frequency$CTRL_0
snv_frequency$GBH_CTRL_T2 <- snv_frequency$GBH_28 - snv_frequency$CTRL_28

#snv_frequency$GBH_CTRL_change_abs <- abs(snv_frequency$GBH_CTRL_change)
#snv_frequency$GBH_change_abs <- abs(snv_frequency$GBH_change)
#snv_frequency$CTRL_change_abs <- abs(snv_frequency$CTRL_change)
snv_frequency$GBH_CTRL_T1_abs <- abs(snv_frequency$GBH_CTRL_T1)
#snv_frequency$GBH_CTRL_T2_abs <- abs(snv_frequency$GBH_CTRL_T2)

snv_frequency$shift <- with(snv_frequency, ifelse((GBH_CTRL_T1_abs <= 0.3 & GBH_change <= -0.7 & GBH_CTRL_T2 <= -0.7), "sig_shift", "not_sig"))
write.csv(snv_frequency, "data files/snv_frequency_changes_07.csv", row.names = F)

## get gene info

mag_list <- c("MAG_00097_1", "MAG_00110_1", "MAG_00179_1", "MAG_00194_1", "MAG_00197_1", "MAG_00201_1", "MAG_00674_1")

sample_names <- read_csv("data files/chapter_1_sample_names.csv")
mag_scaffold_info <- read_tsv("data files/T1_refined.stb", col_names = F)
colnames(mag_scaffold_info) <-  c("scaffold", "mag")

gene_files <- list.files("data files/T1_refined_inStrain", recursive = T, pattern = ".*gene_info.tsv", full.names = T)

all_genes <- data.frame()
for(i in 1:length(gene_files)){
  pond_time_genes <- read.table(gene_files[i], sep="\t", header = T)
  pond_time_genes$Sample <- gene_files[i] %>% substr(32,43)
  all_genes <- rbind(all_genes, pond_time_genes)
}

all_genes <- left_join(all_genes, sample_names)
all_genes <- left_join(all_genes, mag_scaffold_info)
all_genes$mag <- all_genes$mag %>% str_sub(end = -4)
all_genes$Time <- ifelse(grepl("_1", all_genes$Pond), "Day 0", "Day 28")
all_genes$Name_Time <- paste(all_genes$Time, all_genes$Name, sep = " ")
all_genes <- all_genes[, c(24,1,2,26,3:5,8:10,14,17)]
all_genes <- all_genes %>% rename("gene_coverage" = "coverage")
all_genes <- all_genes %>% rename("gene_breadth" = "breadth")
#write.csv(all_genes, "data files/T1_MAG_genes.csv", row.names = F)

mag_genes <- subset(all_genes, mag %in% mag_list)
mag_info <- read_csv("data files/T1_mag_info.csv")
mag_info$time <- with(mag_info, ifelse(time == "1", "Day 0", "Day 28"))
mag_info$Name_Time <- paste(mag_info$time, mag_info$Name, sep = " ")
mag_info <- mag_info[, c(1:3,10)]
mag_genes <- left_join(mag_genes, mag_info)
mag_genes <- subset(mag_genes, mag_coverage >= 5 & mag_breadth >= 0.7)
#write.csv(mag_genes, "data files/MAG_gene_info.csv", row.names = F)

mag_genes$rel_cov <- mag_genes$gene_coverage / mag_genes$mag_coverage
mag_genes <- subset(mag_genes, rel_cov <= 1.2 & rel_cov >= 0.6)

sig_changes <- subset(snv_frequency, shift == "sig_shift")
sig_changes <- subset(sig_changes, gene %in% mag_genes$gene)
sig_nonsyn_changes <- subset(sig_changes, mutation_type == "N")
sig_syn_changes <- subset(sig_changes, mutation_type == "S")

sig_snvs_sum <- sig_changes %>% group_by(mag, gene) %>% count()
sig_snvs_sum <- subset(sig_snvs_sum, !str_detect(gene, ","))

sig_nonsyn_sum <- sig_nonsyn_changes %>% group_by(mag, gene) %>% count()
sig_syn_sum <- sig_syn_changes %>% group_by(mag, gene) %>% count()

eggnog_genes <- read_tsv("data files/eggnog_genes.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("gene", "COG_ID", "Description", "Preferred_name")]
eggnog_genes$scaffold <- eggnog_genes$gene %>% str_extract("[^_]*_[^_]*")
background_cog_snv <- subset(eggnog_genes, gene %in% mag_genes$gene)
background_cog_snv <- subset(background_cog_snv, is.na(COG_ID) == F)
#write.csv(background_cog_snv, "data files/cog_background_snv_freq_genes.csv", row.names = F)

significant_genes <- left_join(sig_snvs_sum, background_cog_snv) 
write.csv(significant_genes, "data files/allele_shifts_significant_genes_07.csv", row.names = F)

significant_nonsyn_genes <- left_join(sig_nonsyn_sum, background_cog_snv) 
write.csv(significant_nonsyn_genes, "data files/nonsyn_allele_shifts_significant_genes_07.csv", row.names = F)

significant_syn_genes <- left_join(sig_syn_sum, background_cog_snv) 
write.csv(significant_syn_genes, "data files/syn_allele_shifts_significant_genes_07.csv", row.names = F)

sum_genes <- significant_genes %>% group_by(mag) %>% count()
write.csv(sum_genes, "data files/sum_sig_shifts_07.csv", row.names = F)

sum_syn_genes <- significant_syn_genes %>% group_by(mag) %>% count()
write.csv(sum_syn_genes, "data files/sum_syn_sig_shifts_07.csv", row.names = F)

sum_nonsyn_genes <- significant_nonsyn_genes %>% group_by(mag) %>% count()
write.csv(sum_nonsyn_genes, "data files/sum_nonsyn_sig_shifts_07.csv", row.names = F)

aro_genes <- significant_nonsyn_genes[grepl("aro", significant_nonsyn_genes$Preferred_name),]

                    