library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00110_1", "MAG_00179_1", "MAG_00194_1", "MAG_00197_1", "MAG_00201_1","MAG_00674_1")

all_samples <- read_csv("data files/chapter_1_sample_names.csv")
all_samples$Time <- ifelse(grepl("_1", all_samples$Pond), "Day 0", "Day 28")
all_samples <- subset(all_samples, !(Name == "CTRL E"))

all_snv <- read_csv("data files/MAG_SNV_depth_info.csv")
all_snv$Time <- ifelse(grepl("1", all_snv$time), "Day 0", "Day 28")
all_snv$Treatment <- ifelse(grepl("CTRL", all_snv$Name), "Control", "GBH")
mag_snv <- subset(all_snv, mag_coverage >= 5 & mag_breadth >= 0.7)

all_sum <- read_csv("data files/T1_SNV_summary_MAG.csv")
all_sum <- subset(all_sum, mag %in% mag_list & mag_coverage >= 5 & mag_breadth >= 0.7)

for(i in 1:length(mag_list)){
  bin_sum <- subset(all_sum, mag == mag_list[i])
  bin_sum <- full_join(bin_sum, all_samples[, c(1,3,4)])
  bin_sum$Treatment <- ifelse(grepl("CTRL", bin_sum$Name), "Control", "GBH")
  
  bin_snv <- subset(mag_snv, mag == mag_list[i])
  bin_snv <- full_join(bin_snv, all_samples[, c(1,3,4)])
  bin_snv <- subset(bin_snv, final_ref_freq < 1 & final_ref_freq > 0)
  bin_snv <- bin_snv %>% complete(group, Name, Time)
  bin_snv <- bin_snv %>% group_by(Name) %>% fill(Treatment, .direction = "updown")
  bin_snv$Name <- factor(bin_snv$Name, levels = c("CTRL A", "CTRL B", "CTRL C", "CTRL D", "GBH A", "GBH B", "GBH C", "GBH D"))
  
  bin_snv_sum <- ggplot(bin_sum, aes(x = Name, y = SNVs_Mbp, fill = Treatment))+
    geom_col()+
    theme_classic()+
    scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
    labs(x = "", y = "Polymorphic sites per Mbp", fill = "")+
    theme(text = element_text(size = 8), axis.text = element_text(colour = "black", size = 6), strip.text.x = element_text(size = 8), legend.position = "none")+
    facet_wrap(~Time, scales = "fixed")
  
  bin_sns_sum <- ggplot(bin_sum, aes(x = Name, y = SNSs_Mbp, fill = Treatment))+
    geom_col()+
    theme_classic()+
    scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
    labs(x = "", y = "Fixed subsitutions per Mbp", fill = "Treatment")+
    theme(text = element_text(size = 8), axis.text = element_text(colour = "black", size = 6), strip.text.x = element_blank(), legend.position = "none")+
    facet_wrap(~Time, scales = "fixed")
  
  bin_snvs <- ggplot(bin_snv, aes(x = Name, y = final_ref_freq, colour = Treatment)) +
    geom_boxplot(na.rm = F, outliers = T,outlier.size = 0.5)+
    scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
    ylim(0,1)+
    theme_classic()+
    scale_x_discrete(drop = FALSE)+
    labs(x = "Pond", y = "Reference frequency of SNVs")+
    theme(text = element_text(size = 8), axis.text = element_text(colour = "black", size = 6), strip.text.x = element_blank(), legend.position = "bottom")+
    facet_wrap(~Time, scales = "fixed")
  
  bin_plot <- bin_snv_sum / bin_sns_sum / bin_snvs
  ggsave(paste("figures/", mag_list[i], "_snvs.pdf", sep = ""), bin_plot, units = "cm", width = 17, height = 17)
}
