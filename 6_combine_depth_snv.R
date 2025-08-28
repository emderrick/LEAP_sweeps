library(tidyverse)
library(dplyr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

depth_files <- list.files("data files/refined_mag_depth/", recursive = T, pattern = "depth_info.csv", full.names = T)
mag_scaffold_info <- read_tsv("data files/T1_refined.stb", col_names = F)
colnames(mag_scaffold_info) <-  c("scaffold", "mag")
mag_SNVs <- read_csv("data files/good_refined_MAG_SNVs.csv")

all_depth <- data_frame()
for(i in 1:length(depth_files)){
  sample_depth <- read_csv(depth_files[i])
  sample_depth$Sample <- depth_files[i] %>% substr(31,42)
  sample_depth <- left_join(sample_depth, mag_scaffold_info)
  sample_depth$position <- str_pad(sample_depth$position, 7, pad = "0")
  sample_depth$mag <- sample_depth$mag %>% str_remove(".fa")
  sample_depth$group <- paste(sample_depth$mag, sample_depth$scaffold, sample_depth$position, sep = "_")
  sample_depth <- sample_depth[, c(4,6,3)]
  snv_depth <- subset(sample_depth, group %in% mag_SNVs$group)
  all_depth <- rbind(all_depth, snv_depth)
}

write.csv(all_depth, "data files/T1_mag_depth_info.csv", row.names = F)

mag_SNV_depth <- full_join(mag_SNVs, all_depth)
mag_SNV_depth$new_ref_freq <- with(mag_SNV_depth, ifelse(depth >= 5, 1, NA))
mag_SNV_depth$final_ref_freq <- with(mag_SNV_depth, ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
mag_SNV_depth <- mag_SNV_depth %>% group_by(group) %>% fill(scaffold, gene, mag, mag_length, length, .direction = "updown")
mag_SNV_depth <- mag_SNV_depth %>% group_by(Sample) %>% fill(Name, time, .direction = "updown")
mag_SNV_depth <- mag_SNV_depth %>% group_by(Sample, mag) %>% fill(mag_coverage, mag_breadth, .direction = "updown")
mag_SNV_depth <- mag_SNV_depth %>% ungroup()
mag_SNV_depth <- subset(mag_SNV_depth, !(is.na(mag_coverage)))
write.csv(mag_SNV_depth, "data files/MAG_SNV_depth_info.csv", row.names = F)


