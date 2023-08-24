library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

MAG_sweep_wide <- data.frame()
for(MAG in mag_list){
  MAG_sweep  <- read.csv(paste("sweep files wide/", MAG, "_sweep_wide.csv", sep = ""), check.names = F)
  MAG_sweep_wide <- bind_rows(MAG_sweep_wide, MAG_sweep)
}

write.csv(MAG_sweep_wide, "all_MAG_sweep_wide.csv")

MAG_sweep_wide$mag <- MAG_sweep_wide$groups %>% substr(1,12)
MAG_sweep_wide$position <- MAG_sweep_wide$groups %>% substr(27,33) %>% str_remove("^0+")
MAG_sweep_wide$control_mean <- rowMeans(MAG_sweep_wide[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], na.rm = T)
MAG_sweep_wide$GBH_mean <- rowMeans(MAG_sweep_wide[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], na.rm=T)
MAG_sweep_wide$abs_val <- abs(MAG_sweep_wide$control_mean - MAG_sweep_wide$GBH_mean)
MAG_sweep_wide$abs_val <- with(MAG_sweep_wide, ifelse(is.nan(abs_val), NA, abs_val))

threshold_snvs <- read_csv("threshold_snvs.csv")

MAG_sweep_wide <- left_join(MAG_sweep_wide, threshold_snvs[,c('groups', 'pass')])
MAG_sweep_wide$pass <- with(MAG_sweep_wide, ifelse(is.na(pass), 'no', pass))

scaffold_info <- read.csv("ANI_95_all_scaffolds.csv") #some scaffolds are missing the length (ones that weren't in SNV file)
scaffold_info <- scaffold_info[,c('scaffold', 'length')] %>% distinct()
unique_sweep <- left_join(MAG_sweep_wide[,c('scaffold', 'mag')], scaffold_info, by = "scaffold") %>% group_by(scaffold) %>% fill(length, .direction = "updown") %>% distinct() 
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

MAG_sweep_plot <- left_join(MAG_sweep_wide[,c('mag', 'scaffold', 'position', 'groups', 'abs_val', 'pass')], unique_sweep_new_pos[,c('scaffold', 'scaffold_position')], by = 'scaffold')
MAG_sweep_plot$abs_val_pass <- with(MAG_sweep_plot, ifelse(pass=="yes", abs_val, NA))
write.csv(MAG_sweep_plot, "MAG_sweep_plot.csv")

test <- MAG_sweep_plot %>% group_by(scaffold) %>% colMeans(abs_val, nrow=1000)

L8_MAG_00011_sweep <- ggplot(subset(MAG_sweep_plot, mag == "L8_MAG_00011"), aes(x = reorder(groups, scaffold_position), y = abs_val, group = mag))+
  geom_line()+
  geom_point(data = filter(MAG_sweep, pass == "yes"), colour = "red", size = 3)+
  geom_vline(xintercept = unique(MAG_sweep$scaffold_position), colour = "blue")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme_classic()

ggsave("L8_MAG_00011_sweep.png", limitsize = F, width = 64, height = 12)


