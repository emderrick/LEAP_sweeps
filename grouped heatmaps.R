library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00054_1", "MAG_00084_1", "MAG_00097_1", "MAG_00103_1",
              "MAG_00110_1", "MAG_00179_1", "MAG_00194_1", "MAG_00197_1",
              "MAG_00201_1", "MAG_00247", "MAG_00636_1", "MAG_00674_1")

all_snv <- read_csv("refined data files/refined_mag_SNV_depth_info.csv")
all_snv$time <- ifelse(grepl("1", all_snv$time), "Day 0", "Day 28")
mag_snv_avg <- read_csv("refined data files/refined_all_avg_freq.csv")
all_snv <- left_join(all_snv, mag_snv_avg)
all_snv$Treatment <- ifelse(grepl("CTRL", all_snv$Name), "Control", "GBH")
mag_snv <- subset(all_snv, mag_coverage >= 5)

all_sum <- read_csv("refined data files/T1_refined_SNV_summary_MAG.csv")
all_sum <- subset(all_sum, mag %in% mag_list & mag_coverage >= 5)

for(i in 1:length(mag_list)){
  bin_snv <- subset(mag_snv, mag == mag_list[i])
  bin_sum <- subset(all_sum, mag == mag_list[i])
  ponds <- nrow(bin_sum)
  
  bin_snv_sum <- ggplot(bin_sum, aes(x = Name, y = SNVs_Mbp, fill = Treatment))+
    geom_col()+
    theme_bw()+
    scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
    labs(x = "", y = "Polymorphic sites per Mbp", fill = "")+
    theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_text(size = 18), legend.position = "none")+
    facet_wrap(~Time, scales = "free_x")
  
  bin_sns_sum <- ggplot(bin_sum, aes(x = Name, y = SNSs_Mbp, fill = Treatment))+
    geom_col()+
    theme_bw()+
    scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
    labs(x = "", y = "Fixed subsitutions per Mbp", fill = "Treatment")+
    theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank())+
    facet_wrap(~Time, scales = "free_x")
  
  bin_snvs <- ggplot(subset(bin_snv, final_ref_freq < 1 & final_ref_freq > 0), aes(x = Name, y = final_ref_freq, colour = Treatment)) +
    geom_boxplot(outliers = T,outlier.size = 0.5 )+
    scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
    ylim(0,1)+
    theme_bw()+
    labs(x = "Pond", y = "Reference frequency of SNVs")+
    theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank(), legend.position = "none")+
    facet_wrap(~time, scales = "free_x")
  
  bin_plot <- bin_snv_sum / bin_sns_sum / bin_snvs
  ggsave(paste("refined figures/", mag_list[i], "_snvs.pdf", sep = ""), bin_plot, limitsize = F, width = 2*ponds, height = 12)
}
