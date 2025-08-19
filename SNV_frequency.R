library(tidyverse)
library(dplyr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00103_1", "MAG_00110_1", "MAG_00179_1",
              "MAG_00194_1", "MAG_00197_1", "MAG_00201_1","MAG_00647_1")

all_mag_SNVs <- read_csv("refined data files/refined_mag_SNV_depth_info.csv")
all_mag_SNVs <- subset(all_mag_SNVs, mag %in% mag_list)

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

write.csv(mag_snvs_wide, "refined data files/refined_all_frequency.csv", row.names = F)
