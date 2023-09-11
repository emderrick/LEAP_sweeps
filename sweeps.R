library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

MAG_sweep_wide <- tibble()

for(MAG in mag_list){
  MAG_sweep  <- read_csv(paste("sweep files wide/", MAG, "_sweep_wide.csv", sep = ""))
  MAG_sweep_wide <- bind_rows(MAG_sweep_wide, MAG_sweep)
}

MAG_sweep_wide$mag <- MAG_sweep_wide$groups %>% substr(1,12)
MAG_sweep_wide$position <- MAG_sweep_wide$groups %>% substr(27,33) %>% str_remove("^0+")
MAG_sweep_wide$control_mean <- rowMeans(MAG_sweep_wide[, c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], na.rm = T)
MAG_sweep_wide$GBH_mean <- rowMeans(MAG_sweep_wide[, c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], na.rm = T)
MAG_sweep_wide$abs_val <- abs(MAG_sweep_wide$control_mean - MAG_sweep_wide$GBH_mean)
MAG_sweep_wide$abs_val <- with(MAG_sweep_wide, ifelse(is.nan(abs_val), NA, abs_val))
write.csv(MAG_sweep_wide, "MAG_sweep_wide.csv", row.names = F)

threshold_snvs <- read_csv("threshold_snvs.csv")

MAG_sweep_wide <- left_join(MAG_sweep_wide, threshold_snvs[, c('groups', 'pass')])
MAG_sweep_wide$pass <- with(MAG_sweep_wide, ifelse(is.na(pass), 'no', pass))

scaffold_info <- read.csv("ANI_95_all_scaffolds.csv") 
scaffold_info <- scaffold_info[,c('scaffold', 'length')] %>% distinct()
unique_sweep <- left_join(MAG_sweep_wide[, c('scaffold', 'mag')], scaffold_info, by = "scaffold") %>% group_by(scaffold) %>% fill(length, .direction = "updown") %>% distinct() 
#there is a scaffold missing in this (cov must have been low) so I had to manually get the length from the fasta file of the MAG with "cat L3_MAG_00058.fa | grep -A1 "L3_MAG_00058_000000000198" | grep -v ">L3_MAG_00058_000000000198" | wc -c" = 15931
unique_sweep$length[unique_sweep$scaffold == "L3_MAG_00058_000000000198"] <- 15931
unique_sweep <- unique_sweep %>% group_by(mag) %>% arrange(desc(length), .by_group = T)

unique_sweep_new_pos <- data.frame()
for(MAG in mag_list){
  scaffold_position = 0
  unique_MAG <- subset(unique_sweep, mag == MAG)
  for(contig in 1: nrow(unique_MAG)){
    scaffold_position = scaffold_position + unique_MAG$length[contig]
    unique_MAG$scaffold_position[contig] <- scaffold_position
  }
  unique_sweep_new_pos <- bind_rows(unique_sweep_new_pos, unique_MAG)
}

MAG_sweep_plot <- left_join(MAG_sweep_wide[, c('mag', 'scaffold', 'position', 'groups', 'abs_val', 'pass')], unique_sweep_new_pos[, c('scaffold', 'scaffold_position')], by = 'scaffold')
MAG_sweep_plot$abs_val_pass <- with(MAG_sweep_plot, ifelse(pass == "yes", abs_val, NA))
write.csv(MAG_sweep_plot, "MAG_sweep_pass_plot.csv", row.names = F)

for(MAG in mag_list){
  MAG_sweep <- subset(MAG_sweep_plot, mag == MAG)
  sweep_plot <- ggplot(MAG_sweep, aes(x = group, y = abs_val_pass, group = 1))+
  geom_line()+
  geom_vline(xintercept = MAG_sweep$scaffold_position, colour = "blue")+
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(paste(MAG, "_big_sweep.png", sep = ""), limitsize = F, width = 150, height = 12)
}


#to do by bins of 1000
#MAG_sweep_wide$bin <- ((as.numeric(MAG_sweep_wide$position)) - 1) /1000
#MAG_sweep_wide$bin <- floor(MAG_sweep_wide$bin)
#MAG_sweep_wide$one <- 1
#MAG_sweep_sum <- MAG_sweep_wide %>% group_by(mag, scaffold, bin) %>% summarize(sum_abs_val = sum(abs_val, na.rm = T), bin_count = sum(one))
#MAG_sweep_sum$abs_val_mean <- MAG_sweep_sum$sum_abs_val / MAG_sweep_sum$bin_count
#MAG_sweep_sum$one <- 1
#MAG_sweep_scaf <- MAG_sweep_sum %>% group_by(mag, scaffold) %>% summarize(bins_per_scaf = sum(one))
#MAG_sweep_scaf <- MAG_sweep_scaf %>% arrange(desc(bins_per_scaf), .by_group = T)
#sweep_new_pos <- data.frame()
#for(MAG in mag_list){
  #scaffold_position = 0
  #MAG_sweep <- subset(MAG_sweep_scaf, mag == MAG)
  #for(contig in 1: nrow(MAG_sweep)){
    #scaffold_position = scaffold_position + MAG_sweep$bins_per_scaf[contig]
    #MAG_sweep$scaffold_position[contig] <- scaffold_position
  #}
  #sweep_new_pos <- bind_rows(sweep_new_pos, MAG_sweep)
#}
#sweep_new_pos <- sweep_new_pos %>% group_by(mag) %>% mutate(index = row_number())
#MAG_sweep_plot <- left_join(MAG_sweep_sum[, c('mag', 'scaffold', 'abs_val_mean', 'bin', 'bin_count')], sweep_new_pos[, c('scaffold', 'scaffold_position', 'bins_per_scaf', 'index')], by = 'scaffold')
#MAG_sweep_plot$group <- paste(MAG_sweep_plot$bin, MAG_sweep_plot$bins_per_scaf)
#MAG_sweep_plot <- MAG_sweep_plot %>% arrange(index)
#write.csv(MAG_sweep_plot, "MAG_sweep_plot.csv")
#for(MAG in mag_list){
    #MAG_sweep <- subset(MAG_sweep_plot, mag == MAG)
    #sweep_plot <- ggplot(MAG_sweep, aes(x = group, y = abs_val_mean, group = 1))+
    #geom_line()+
    #geom_vline(xintercept = MAG_sweep$scaffold_position, colour = "blue")+
    #theme_classic()+
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    #ggsave(paste(MAG, "_big_sweep.png", sep = ""), limitsize = F, width = 150, height = 12)
#}
                       