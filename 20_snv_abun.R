library(tidyverse)
library(rstatix)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

rel_abun <- read_csv("data files/MAG_rel_abun_change.csv")
mag_snvs <- read_csv("data files/T1_subsamp_SNV_summary_MAG.csv")

SNV_files <- list.files("data files/", recursive = T, pattern = "T1_filt_SNVs", full.names = T)

all_MAG_SNVs <- data_frame()
 for(i in 1:length(SNV_files)){
   all_SNVs <- read_csv(SNV_files[i])
   all_SNVs <- subset(all_SNVs, mag_coverage >= 5, mag_breadth >= 0.7)
   all_SNVs <- all_SNVs[, c("mag", "scaffold", "group", "Name", "time", "position_coverage", "mutation_type", "ref_freq")]
   all_MAG_SNVs <- rbind(all_MAG_SNVs, all_SNVs)
 }

all_MAG_SNVs$time <- as.character(all_MAG_SNVs$time)
write.csv(all_MAG_SNVs, "data files/all_MAG_SNVs.csv", row.names = F)

mag_snvs$total_variants <- mag_snvs$SNSs_Mbp + mag_snvs$SNVs_Mbp
mag_snvs$N_total_variants <- mag_snvs$N_SNSs_Mbp + mag_snvs$N_SNVs_Mbp

mag_snv_change <- mag_snvs[, c(1,21,17:19)]
mag_snv_change <- pivot_wider(mag_snv_change, names_from = Time, values_from = total_variants)
mag_snv_change <- mag_snv_change[complete.cases(mag_snv_change ), ]
mag_snv_change$snv_change <- mag_snv_change$`Day 28` - mag_snv_change$`Day 0`
mag_snv_change$mag_pond <- paste(mag_snv_change$mag, mag_snv_change$Name)

mag_n_snv_change <- mag_snvs[, c(1,22,17:19)]
mag_n_snv_change <- pivot_wider(mag_n_snv_change, names_from = Time, values_from = N_total_variants)
mag_n_snv_change <- mag_n_snv_change[complete.cases(mag_n_snv_change ), ]
mag_n_snv_change$N_snv_change <- mag_n_snv_change$`Day 28` - mag_n_snv_change$`Day 0`

abun_snv_all <- left_join(rel_abun[, c(1:3,6:12)], mag_snv_change[,c(1,3,6)])
abun_snv_all <- left_join(abun_snv_all, mag_n_snv_change[, c(1,3,6)])
abun_snv_all <- abun_snv_all[complete.cases(abun_snv_all), ]
abun_snv_all$arg <- with(abun_snv_all, ifelse(total_hits >= 1, "Yes", "No"))
abun_snv_all$mag_pond <- paste(abun_snv_all$mag, abun_snv_all$Name)

all_MAG_SNVs$mag_pond <- paste(all_MAG_SNVs$mag, all_MAG_SNVs$Name)
all_MAG_SNVs <- subset(all_MAG_SNVs, mag_pond %in% abun_snv_all$mag_pond)
write.csv(all_MAG_SNVs, "data files/snvs_for_depth.csv", row.names = F)

mag_scaffolds <- as.data.frame(all_MAG_SNVs$scaffold)
mag_scaffolds <- unique(mag_scaffolds)
colnames(mag_scaffolds) <- "scaffold"
write_tsv(mag_scaffolds, "data files/all_mag_scaffolds.txt")


depth_files <- list.files("data files/inStrain all mag depth", recursive = T, pattern = "depth.csv", full.names = T)
mag_scaffold_info <- read_tsv("data files/T1_refined.stb", col_names = F)
colnames(mag_scaffold_info) <-  c("scaffold", "mag")

for(i in 1:length(depth_files)){
  sample_depth <- read_csv(depth_files[i])
  sample_depth <- left_join(sample_depth, mag_scaffold_info)
  sample_depth$mag <- sample_depth$mag %>% str_remove(".fa")
  sample_depth$position <- str_pad(sample_depth$position, 7, pad = "0")
  sample_depth$group <- paste0(sample_depth$mag, "_", sample_depth$scaffold, "_", sample_depth$position)
  sample_depth$Sample <- depth_files[i] %>% substr(35,46)
  sample_depth <- sample_depth[, c(6,5,1)]
  snv_depth <- subset(sample_depth, group %in% all_MAG_SNVs$group)
  sample <- depth_files[i] %>% substr(35,46)
  write.csv(snv_depth, paste("data files/MAG_SNV_depth_info_", sample, ".csv", sep = ""), row.names = F)
}



