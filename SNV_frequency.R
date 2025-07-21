library(tidyverse)
library(dplyr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("bin.305", "bin.384", "bin.609", "bin.626", "bin.651", "bin.663",
              "bin.676", "bin.707", "bin.870")

all_mag_SNVs <- read_csv("data files/mag_SNV_depth_info.csv")
all_mag_SNVs$time <- with(all_mag_SNVs, ifelse(time == "1", "Day 0", "Day 28"))
all_mag_SNVs$pond_time <- paste(all_mag_SNVs$time, all_mag_SNVs$Name, sep = " ")

all_avg_freq <- data_frame()
all_frequency <- data_frame()
for(i in 1:length(mag_list)){
  mag_snvs <- subset(all_mag_SNVs, mag == mag_list[i])
  mag_snvs <- mag_snvs[, c("gene", "pond_time", "mag", "group", "final_ref_freq")]
  mag_snvs_wide <- pivot_wider(mag_snvs, names_from = pond_time, values_from = final_ref_freq)
  mag_snvs_wide$sum_na <- rowSums(is.na(mag_snvs_wide))
  mag_snvs_wide$CTRL_0 <-rowMeans(mag_snvs_wide[,grep("Day 0 CTRL",colnames(mag_snvs_wide))], na.rm = T)
  mag_snvs_wide$GBH_0 <- rowMeans(mag_snvs_wide[,grep("Day 0 GBH",colnames(mag_snvs_wide))], na.rm = T)
  mag_snvs_wide$CTRL_28 <- rowMeans(mag_snvs_wide[,grep("Day 28 CTRL",colnames(mag_snvs_wide))], na.rm = T)
  mag_snvs_wide$GBH_28 <- rowMeans(mag_snvs_wide[,grep("Day 28 GBH",colnames(mag_snvs_wide))], na.rm = T)
  mag_snvs_wide <- mag_snvs_wide %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
  mag_snvs_wide$med <- rowMeans(mag_snvs_wide[,grep("Day",colnames(mag_snvs_wide))], na.rm = T)
  mag_snvs_wide$missing <- rowSums(is.na(mag_snvs_wide[, c("CTRL_0", "GBH_0", "CTRL_28", "GBH_28")]))
  mag_snvs_avg <- pivot_longer(mag_snvs_wide, cols = grep("Day", colnames(mag_snvs_wide)), names_to = "Name_time", values_to = "Reference_Frequency" )
  mag_snvs_avg <- mag_snvs_avg[, c("mag", "group", "med")] %>% unique()
  all_avg_freq <- rbind(all_avg_freq, mag_snvs_avg)
  mag_snvs_wide <- mag_snvs_wide[, c("mag", "gene", "group", "CTRL_0", "GBH_0", "CTRL_28", "GBH_28", "med")]
  all_frequency <- rbind(all_frequency, mag_snvs_wide)
}

write.csv(all_frequency, "data files/all_frequency.csv", row.names = F)
write.csv(all_avg_freq, "data files/all_avg_freq.csv", row.names = F)
