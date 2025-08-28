library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00110_1", "MAG_00179_1", "MAG_00194_1",
              "MAG_00197_1", "MAG_00201_1", "MAG_00674_1")

SNV_files <- list.files("data files/", recursive = T, pattern = "filt_SNVs", full.names = T)

mag_SNVs <- data_frame()
for(i in 1:length(SNV_files)){
  all_SNVs <- read_csv(SNV_files[i])
  mag_only <- subset(all_SNVs, mag %in% mag_list)
  mag_SNVs <- rbind(mag_SNVs, mag_only)
}

write.csv(mag_SNVs, "data files/good_refined_MAG_SNVs.csv", row.names = F)

mag_scaffolds <- as.data.frame(mag_SNVs$scaffold)
mag_scaffolds <- unique(mag_scaffolds)
colnames(mag_scaffolds) <- "scaffold"
write.csv(mag_scaffolds, "data files/good_refined_mag_scaffolds.csv", row.names = F)
