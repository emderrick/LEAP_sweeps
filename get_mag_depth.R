library(tidyverse)

depth_files <- list.files("/mfs/ederrick/chapter_1/08_T1_T3/", recursive = T, pattern = "depth.txt", full.names = T)
scaffold_list <- read_csv("good_MAG_scaffolds.csv")

for(i in 1:length(depth_files)){
  sample <- depth_files[i] %>% substr(35,46)
  all_depth <- read_table(depth_files[i])
  colnames(all_depth) <- c("scaffold", "position", "depth")
  mag_depth <- subset(all_depth, scaffold %in% scaffold_list$scaffold)
  write.csv(mag_depth, paste(sample, "_depth_info.csv", sep = ""), row.names = F)
}