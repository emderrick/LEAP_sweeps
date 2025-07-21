library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("bin.259", "bin.305", "bin.384", "bin.609", "bin.626", "bin.651", "bin.663",
              "bin.676", "bin.707", "bin.870")

SNV_files <- list.files("data files/", recursive = T, pattern = "T1_filtered_SNVs", full.names = T)

mag_SNVs <- data_frame()
for(i in 1:length(SNV_files)){
  sample <- SNV_files[i] %>% substr(30,41)
  all_SNVs <- read_csv(SNV_files[i])
  mag_only <- subset(all_SNVs, mag %in% mag_list)
  mag_SNVs <- rbind(mag_SNVs, mag_only)
}

mag_SNVs <- subset(mag_SNVs, mag_coverage >= 5 & mag_breadth >= 0.5)
write.csv(mag_SNVs, "data files/good_MAG_SNVs.csv", row.names = F)

mag_scaffolds <- as.data.frame(mag_SNVs$scaffold)
mag_scaffolds <- unique(mag_scaffolds)
colnames(mag_scaffolds) <- "scaffold"
write.csv(mag_scaffolds, "data files/good_mag_scaffolds.csv", row.names = F)
