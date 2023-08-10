library(tidyverse)
library(dplyr)

filtered_SNVs <- read_csv("filtered_ANI_95_mag_SNVs.csv")

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

ponds_list <- list(I4_MAG_00006 = list("Control B at T2", "Control E at T2", "GBH A at T2", "GBH D at T2"), 
                   I4_MAG_00065 = list("Control A at T2", "Control B at T2", "Control C at T2", "Control D at T2", "Control E at T2", "GBH B at T2"), 
                   L2_MAG_00052 = list("Control A at T2", "Control B at T2", "Control D at T2", "Control E at T2", "GBH A at T2"), 
                   L3_MAG_00058 = list("Control C at T2", "Control D at T2", "GBH C at T2", "GBH D at T2"), 
                   L4_MAG_00099 = list("Control D at T2", "GBH A at T2", "GBH B at T2", "GBH C at T2", "GBH D at T2"), 
                   L7_MAG_00020 = list("Control A at T1", "Control A at T2", "Control C at T1", "Control C at T2", "Control D at T1", "Control D at T2", "GBH A at T1", "GBH A at T2"), 
                   L7_MAG_00028 = list("Control E at T2", "GBH A at T2", "GBH B at T2", "GBH C at T2"), 
                   L7_MAG_00043 = list("Control D at T2", "GBH B at T2", "GBH C at T2"), 
                   L8_MAG_00011 = list("Control E at T2", "GBH A at T2", "GBH D at T2"), 
                   L8_MAG_00019 = list("Control E at T2", "GBH A at T2", "GBH B at T2", "GBH D at T2"), 
                   L8_MAG_00042 = list("Control A at T2", "Control C at T2", "Control D at T2", "Control E at T2", "GBH D at T2"))

merge_depth <- function(filtered_SNVs, MAG, ponds_list){
  MAG_snv <- filter(filtered_SNVs, mag==MAG)
  MAG_snv$groups <- paste(MAG_snv$scaffold, str_pad(MAG_snv$position, 7, pad = "0"))
  MAG_snv <-  MAG_snv[,c(1:17, 19:22, 41, 45, 72, 74:76, 78:82)]
  full_MAG_snv <- complete(MAG_snv, timepoint, groups)
  
  MAG_depth <- read.table(paste(MAG, "_depth.txt", sep=""), sep="\t", header=F)
  MAG_depth <- MAG_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
  MAG_depth$groups <- paste(MAG_depth$scaffold, str_pad(MAG_depth$position, 7, pad = "0"))
  MAG_depth <- left_join(MAG_depth, MAG_snv[,c(19,27,32)], by= c("timepoint", "groups"))
  MAG_depth <- MAG_depth %>% group_by(timepoint) %>% fill(new_name, .direction = "updown")
  MAG_index<-which(names(ponds_list)==MAG)
  MAG_depth <- subset(MAG_depth, new_name %in% ponds_list[[MAG_index]])
  MAG_depth$new_ref_freq <- with(MAG_depth, ifelse(samtools_depth >= 5, 1, NA))

  all_MAG <- left_join(MAG_depth, MAG_snv, by=c("timepoint", "groups", "scaffold", "position", "new_name"))
  all_MAG$final_ref_freq <- with(all_MAG, ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
  all_MAG <- all_MAG %>% group_by(groups) %>% fill(scaffold, position, gene, mag, mag_length, length, .direction="updown")
  all_MAG <- all_MAG %>% group_by(timepoint) %>% fill(new_time, treatment, name, new_name, full_group, .direction="updown")
  all_MAG <- all_MAG[,c(1:6, 8:12, 35, 14:34)]
  #save files for sweep plots
  write.csv(all_MAG, paste(MAG, "_sweep_snvs.csv", sep=""), row.names=F)
  small_sweep <- all_MAG[,c('scaffold', 'position', 'length', 'groups', 'gene', 'mag', 'mag_length', 'new_name', 'final_ref_freq')]
  write.csv(small_sweep, paste("small_", MAG, "_sweep.csv", sep=""), row.names=F)
  MAG_sweep_wide <- all_MAG[,c(1, 23, 19, 5, 6, 12)] %>% pivot_wider(names_from = new_name, values_from = final_ref_freq)
  write.csv(MAG_sweep_wide, paste(MAG, "_sweep_wide.csv", sep=""), row.names=F)
  print(paste("done", MAG, "sweep files"))
  #save heatmap file
  MAG_heat <- subset(all_MAG, groups %in% full_MAG_snv$groups)
  write.csv(MAG_heat, paste("all_", MAG, "_SNVs.csv", sep=""), row.names=F)
  print(paste("done", MAG, "heatmap file"))
}

for(MAG in mag_list){
  merge_depth(filtered_SNVs, MAG, ponds_list)
}

all_MAG_SNVs <- data.frame()
for(MAG in mag_list){
  MAG_SNVs  <- read_csv(paste("all_", MAG, "_SNVs.csv", sep=""))
  all_MAG_SNVs <- rbind(all_MAG_SNVs, MAG_SNVs)
  print(paste("done", MAG))
}

write.csv(all_MAG_SNVs, "all_MAG_SNVs.csv", row.names=F)
