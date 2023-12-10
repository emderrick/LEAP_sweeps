library(tidyverse)
library(dplyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

filtered_SNVs <- read_csv("filtered_mag_SNVs_subsamp.csv")

depth_files <- list.files("instrain depth", recursive = T, pattern = ".*.csv", full.names = T)
all_depth <- data.frame()

for(i in 1:length(depth_files)){
  pond_time_mags <- read_csv(depth_files[i])
  timepoint <- gsub(".*instrain depth", "", depth_files[i]) %>% substr(2,10)
  pond_time_mags <- cbind(pond_time_mags,timepoint = rep(timepoint, nrow(pond_time_mags)))
  all_depth <- rbind(all_depth, pond_time_mags)
}

all_depth$mag<- all_depth$scaffold %>% substr(1,12)
write.csv(all_depth, "all_instrain_depth.csv", row.names = F)

for(MAG in mag_list){
  MAG_snv <- subset(filtered_SNVs, mag == MAG)
  MAG_snv$groups <- paste(MAG_snv$scaffold, str_pad(MAG_snv$position, 7, pad = "0"))
  MAG_snv <-  MAG_snv[, c(1:17, 19:22, 41, 45, 72, 74:76, 78:82)]
  
  MAG_depth <- subset(all_depth, mag == MAG)
  MAG_depth$groups <- paste(MAG_depth$scaffold, str_pad(MAG_depth$position, 7, pad = "0"))
  MAG_depth <- subset(MAG_depth, groups %in% MAG_snv$groups)
  MAG_depth <- left_join(MAG_depth, MAG_snv[, c(19,27,32)], by = c("timepoint", "groups"))
  MAG_depth <- MAG_depth %>% group_by(timepoint) %>% fill(new_name, .direction = "updown")
  MAG_depth$new_ref_freq <- with(MAG_depth, ifelse(coverage >= 5, 1, NA))
  
  all_MAG <- left_join(MAG_depth, MAG_snv, by = c("timepoint", "groups", "mag", "scaffold", "position", "new_name"))
  all_MAG$final_ref_freq <- with(all_MAG, ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
  all_MAG <- all_MAG[, c(3:7, 9:13, 35, 14:34)]
  all_MAG <- ungroup(all_MAG)
  all_MAG <- complete(all_MAG, timepoint, groups)
  all_MAG <- all_MAG %>% group_by(groups) %>% fill(scaffold, gene, mag, mag_length, length, .direction = "updown")
  all_MAG <- all_MAG %>% group_by(timepoint) %>% fill(new_time, treatment, name, new_name, full_group, .direction = "updown")
  write.csv(all_MAG, paste("all_", MAG, "_SNVs_subsamp.csv", sep = ""), row.names = F)
}

all_MAG_SNVs <- data.frame()
for(MAG in mag_list){
  MAG_SNVs <- read_csv(paste("all_", MAG, "_SNVs_subsamp.csv", sep = ""))
  all_MAG_SNVs <- rbind(all_MAG_SNVs, MAG_SNVs)
}

write.csv(all_MAG_SNVs, "all_MAG_SNVs_subsamp.csv", row.names = F)
