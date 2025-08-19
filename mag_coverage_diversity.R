library(tidyverse)
library(ggplot2)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")
mag_list <- c("bin.305", "bin.609", "bin.676")

mag_pres <- read_csv("data files/T1_SNV_summary_MAG.csv")
mag_pres <- subset(mag_pres, med_cov >= 5 & mag_breadth >= 0.5)
mag_pres$Treatment_Time <- paste(mag_pres$Time, mag_pres$Treatment, sep = " ")

mag_pres <- subset(mag_pres, mag %in% mag_list)

mag_coverage <- ggplot(mag_pres, aes(x = mag_coverage, y = SNVs_Mbp, colour = Treatment))+
  geom_point(size = 2)+
  theme_bw()+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  facet_wrap(~mag, scales = "free")
ggsave("figures/T1_mag_coverage_vs_diversity.pdf", mag_coverage,  limitsize = F, width = 10, height = 3)

mag_diversity_pond <- ggplot(mag_pres, aes(x = Treatment_Time, y = SNVs_Mbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 1, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "SNVs per Mbp")+  
  theme_bw()+
  facet_wrap(~mag, scales = "free")
ggsave("figures/T1_mag_diversity_pond.pdf", mag_diversity_pond, limitsize = F, width = 14, height = 4)


