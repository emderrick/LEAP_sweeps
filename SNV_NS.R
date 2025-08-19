library(tidyverse)
library(ggplot2)
library(rstatix)

options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")
snv_count <- read_csv("data files/T1_limited_SNV_info.csv")

mag_info <- read_csv("data files/T1_mag_info.csv")
mag_SNVs <- read_csv("data files/T1_SNV_summary_MAG.csv")
NS_MAG <- read_csv("data files/T1_NS_summary_MAG.csv")
mag_SNV_NS <- left_join(mag_SNVs, NS_MAG)
mag_SNV_NS_10x <- subset(mag_SNV_NS, med_cov >= 5 & mag_breadth >= 0.4)

SNV_NS_plot <- ggplot(mag_SNV_NS_10x, aes(x = NS_ratio, y = SNVs_Mbp, colour = Time))+
  geom_point(size = 1)+
  scale_colour_manual(values = c("#FDE725FF", "#482677FF"))+
  labs(x = "N:S ratio", y = "Polymorphic sites per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))+
  facet_wrap(~Treatment)

ggsave("figures/SNV_NS_10x.pdf", SNV_NS_plot, limitsize = F, width = 8, height = 3)

