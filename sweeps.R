library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L2_MAG_00052 = "Erythrobacter sp.",
               L3_MAG_00058 = "Prosthecobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B sp.",
               L7_MAG_00028 = "SYFN01 sp.", L7_MAG_00043 = "Luteolibacter sp.",L8_MAG_00011 = "UBA953 sp.",
               L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp."))

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

MAG_sweep_plot <- read_csv("MAG_sweep_pass_plot.csv")


for(MAG in mag_list){
  MAG_sweep <- subset(MAG_sweep_plot, mag == MAG)
  sweep_plot <- ggplot(MAG_sweep, aes(x = groups, y = abs_val_pass))+
    geom_point()+
    geom_vline(xintercept = unique(MAG_sweep$scaffold_position), colour = "blue")+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
    theme_classic()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    labs(y = "Absolute Value Difference", x = "Scaffold")
  save_plot(paste(MAG, "_pass_sweep.png", sep = ""), sweep_plot, base_asp = 6)

}

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs$snv_count <- with(all_MAG_snvs, ifelse(is.na(number_divergent) == T, 0, number_divergent))
MAG_snvs_scaffold_sum <- all_MAG_snvs %>% group_by(mag, scaffold, length, new_name) %>% summarize(total = sum(snv_count))
MAG_snvs_scaffold_sum$snv_per_kb <- (MAG_snvs_scaffold_sum$total/MAG_snvs_scaffold_sum$length) * 1000

for(MAG in mag_list){
  MAG_SNVs <- subset(MAG_snvs_scaffold_sum, mag == MAG)
  rows <- length(unique(MAG_SNVs$new_name))/2
  scaffold_hist <- ggplot(MAG_SNVs, aes(x = snv_per_kb), group = mag)+
    geom_histogram(bins = 25, fill = "grey", colour = "black", linewidth = 0.5)+
    theme_classic()+
    theme(text = element_text(size = 10), axis.text = element_text(colour = 'black'))+
    labs(y = "Frequency", x = "SNVs per kb of scaffold")+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = c(0, 0))+
    facet_wrap(~new_name, ncol = 1)
  
  save_plot(paste(MAG, "_scaffold_hist.jpeg", sep = ""), ncol = 1, nrow = rows, scaffold_hist)
}


threshold_snvs$snv_count <- with(threshold_snvs, ifelse(pass == "yes", 1, 0))
threshold_snvs_sum <- threshold_snvs %>% group_by(mag, scaffold, length) %>% summarize(total = sum(snv_count))
threshold_snvs_sum$snv_per_kb <- (threshold_snvs_sum$total/threshold_snvs_sum$length) * 1000
threshold_snvs_gene_sum <- threshold_snvs %>% group_by(mag, gene, length) %>% summarize(total = sum(snv_count))
threshold_snvs_gene_sum$snv_per_kb <- (threshold_snvs_gene_sum$total/threshold_snvs_gene_sum$length) * 1000

for(MAG in mag_list){
  MAG_SNVs <- subset(threshold_snvs_sum, mag == MAG)
  scaffold_hist <- ggplot(MAG_SNVs, aes(x = snv_per_kb))+
    geom_histogram(bins = 25, fill = "grey", colour = "black", linewidth = 0.5)+
    theme_classic()+
    theme(text = element_text(size = 10), axis.text = element_text(colour = 'black'))+
    labs(y = "Frequency", x = "SNVs per kb of scaffold")+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    scale_x_continuous(expand = c(0, 0))
  
  save_plot(paste(MAG, "_threshold_scaffold_hist.jpeg", sep = ""), scaffold_hist)
}
