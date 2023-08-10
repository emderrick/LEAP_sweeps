library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(patchwork)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099")
                
#"L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

all_sweep <- data.frame()
for(MAG in mag_list){
  sweep  <- read_csv(paste("small_", MAG, "_sweep.csv", sep=""))
  all_sweep <- bind_rows(all_sweep, sweep)
  print(paste("done", MAG))
}

write.csv(all_sweep, "all_MAG_sweep.csv")

sweep_hor <- spread(all_sweep, key=new_name, value=final_ref_freq)
sweep_hor$control_mean <-rowMeans(sweep_hor[,c()], na.rm=T)
sweep_hor$GBH_mean <-rowMeans(sweep_hor[,c()], na.rm=T)
sweep_hor$abs_val<- abs(sweep_hor$control_mean - sweep_hor$GBH_mean)
sweep_hor$abs_val<- with(sweep_hor, ifelse(is.nan(abs_val), NA, abs_val))

threshold_snvs <- read_csv("threshold_snvs.csv")

MAG_sweep <- left_join(sweep_hor, threshold_snvs[,c(4,28)])
MAG_sweep$pass <- with(MAG_sweep, ifelse(is.na(pass), 'no', pass))

unique_sweep <- MAG_sweep[,c('scaffold', 'length')] %>% group_by(scaffold) %>% fill(length, .direction = "updown") %>% distinct() 
unique_sweep <- unique_sweep[order(unique_sweep$length, decreasing = T),]

scaffold_position = 0
for(contig in 1: nrow(unique_sweep)){
  scaffold_position = scaffold_position + unique_sweep$length[contig]
  unique_sweep$scaffold_position[contig] <- scaffold_position
}

MAG_sweep <- left_join(MAG_sweep, unique_sweep[,c(1,3)], by = 'scaffold')

L8_MAG_00011_sweep <- ggplot(MAG_sweep, aes(x = reorder(groups, scaffold_position), y=abs_val, group=mag))+
  geom_line()+
  geom_point(data=filter(MAG_sweep, pass=="yes"), colour="red", size=5)+
  geom_vline(xintercept = unique(MAG_sweep$scaffold_position), colour="blue")+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme_classic2()

ggsave("L8_MAG_00011_sweep.png", limitsize=F, width=64, height=12)




test <- subset(MAG_sweep, scaffold == "L4_MAG_00099_000000000001" | scaffold == "L4_MAG_00099_000000000002" | scaffold == "L4_MAG_00099_000000000003") %>% fill(mag, .direction = "updown") 
test <- test[,c(1:16)]
unique_test <- test[,c('scaffold', 'length')]  %>% group_by(scaffold) %>% fill(length, .direction = "updown") %>% distinct() 
unique_test <- unique_test[order(unique_test$length, decreasing = T),]

scaffold_position = 0
for(contig in 1: nrow(unique_test)){
  scaffold_position = scaffold_position + unique_test$length[contig]
  unique_test$scaffold_position[contig] <- scaffold_position
}

test_sweep <- left_join(test, unique_test[,c(1,3)], by = 'scaffold')


unique_sweep$length[unique_sweep$scaffold == "L4_MAG_00099_000000000214"] <- 3661 #this was missing because there was not an SNV in this scaffold -- got value from ANI_95_all_scaffolds.csv

for(MAG in mag_list){
  sweep  <- read_csv(paste("small_", MAG, "_sweep.csv", sep=""))
  sweep_hor <- spread(sweep, key=new_name, value=final_ref_freq)
  write.csv(sweep_hor, paste(MAG, "_sweep_hor.csv", sep=""))
}
